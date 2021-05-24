// Copyright ©2020 Dan Kortschak. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// smeargol distributes count data across the Gene Ontology DAG provided and
// prints the GO terms, their roots and depths and distributed counts in a
// tsv table to stdout. It logs gene identifiers that do not have GO term
// annotations to stderr. The graph analysis assumes Ensembl gene identifiers
// and Gene Ontology graph structure.
//
// The input counts file is a tab-delimited file with the first column being
// Ensembl gene ID (ENSG00000000000) and remaining columns being count data.
// The first row is expected to be labelled with the first column being Geneid
// and the remaining columns holding the names of the samples.
//
// The Gene Ontology is required to be in Owl format. The file can be
// obtained from http://current.geneontology.org/ontology/go.owl.
//
// The ENSG to GO mapping is expected to be in RDF N-Triples or N-Quads in
// the form:
//
//  <obo:GO_0000000> <local:annotates> <ensembl:ENSG00000000000> .
//
// for each GO term to Ensembl gene annotation.
//
// All input files are expected to be gzip compressed and the output is
// written uncompressed to standard output and standard error.
package main

import (
	"compress/gzip"
	"encoding/csv"
	"flag"
	"fmt"
	"image/color"
	"io"
	"log"
	"math"
	"math/big"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/formats/rdf"
	"gonum.org/v1/gonum/graph/traverse"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"

	"github.com/kortschak/gogo"

	"github.com/kortschak/smeargol/internal/owl"
)

func main() {
	var (
		in       = flag.String("in", "", "specify the counts input (.tsv.gz - required)")
		ontopath = flag.String("ontology", "", "specify the GO file (.owl.gz - required)")
		mappath  = flag.String("map", "", "specify the ENSG to GO mapping (.nt.gz/.nq.gz - required)")
		lean     = flag.Bool("lean", true, "only load relevant parts of ontology")
		cut      = flag.Float64("cut", 1, "minimum valid singular value")
		frac     = flag.Float64("frac", 0.75, "include singular values up to this cumulative fraction")
		debug    = flag.Bool("debug", false, "output binary assignments - only small sets")
		help     = flag.Bool("help", false, "print help text")
	)
	flag.Parse()

	if *help {
		flag.Usage()
		fmt.Fprintf(os.Stderr, `
%s distributes count data across the gene ontology DAG provided and
prints the GO terms, their roots and depths and distributed counts in a
tsv table to stdout. It logs gene identifiers that do not have GO term
annotations to stderr. The graph analysis assumes Ensembl gene identifiers
and Gene Ontology graph structure.

The input counts file is a tab-delimited file with the first column being
Ensembl gene ID (ENSG00000000000) and remaining columns being count data.
The first row is expected to be labelled with the first column being Geneid
and the remaining columns holding the names of the samples.

The Gene Ontology is required to be in Owl format. The file can be
obtained from http://current.geneontology.org/ontology/go.owl.

The ENSG to GO mapping is expected to be in RDF N-Triples or N-Quads in
the form:

 <obo:GO_0000000> <local:annotates> <ensembl:ENSG00000000000> .

for each GO term to Ensembl gene annotation.

All input files are expected to be gzip compressed and the output is
written uncompressed to standard output and standard error.

Copyright ©2020 Dan Kortschak. All rights reserved.

`, filepath.Base(os.Args[0]))
		os.Exit(0)
	}

	if *in == "" || *ontopath == "" || *mappath == "" {
		flag.Usage()
		os.Exit(2)
	}

	log.Println(os.Args)
	for _, d := range []string{
		"matrices",
		"plots",
	} {
		err := os.Mkdir(d, 0o755)
		if err != nil {
			log.Fatal(err)
		}
	}

	log.Println("[loading count data]")
	data, err := mappingCounts(*in)
	if err != nil {
		log.Fatalf("failed to load count data: %v", err)
	}

	if *lean {
		log.Println("[loading lean ontology]")
	} else {
		log.Println("[loading ontology]")
	}
	ontology, err := ontologyGraph(*ontopath, *lean)
	if err != nil {
		log.Fatalf("failed to load ontology: %v", err)
	}

	log.Println("[loading gene to ontology mappings]")
	err = connectGeneIDsTo(ontology, *mappath, data.counts)
	if err != nil {
		log.Fatalf("failed to connect gene IDs to ontology: %v", err)
	}

	log.Println("[smearing counts]")
	roots := ontology.Roots(false)
	sort.Slice(roots, func(i, j int) bool { return roots[i].Value < roots[j].Value })
	ontoData := distributeCounts(ontology, roots, data)

	log.Println("[writing smeared count matrices]")

	var dw *debugWriter
	if *debug {
		dw = newDebugWriter(ontology, ontoData, data)
	}

	var wg sync.WaitGroup
	for i := range ontoData {
		i := i
		wg.Add(1)
		go func() {
			defer wg.Done()
			lastD := -1
			var goTerms []string
			walkDownSubClassesFrom(roots[i], ontology, func(r, t rdf.Term, d int) {
				dw.record(i, d, r, t)

				if lastD == -1 || d == lastD {
					goTerms = append(goTerms, t.Value)
					lastD = d
					return
				}
				lastD = d

				// Write out matrices for this depth. Note that d is now
				// referring to the next level.
				writeCountData(r.Value, d-1, goTerms, data, ontoData[i], *cut, *frac)
				goTerms = goTerms[:0]
				goTerms = append(goTerms, t.Value)
			})

			// Write out last depth.
			writeCountData(roots[i].Value, lastD, goTerms, data, ontoData[i], *cut, *frac)
		}()
	}
	wg.Wait()
	dw.flush()
}

// writeCountData writes out a matrix of gene expression data summed according
// to the bit vector data collected during the walk of the GO DAG. It also
// performs an SVD of the matrix, plotting the singular values and obtaining
// an optimal truncation for each GO level/aspect.
func writeCountData(root string, depth int, goTerms []string, data *countData, ontoData map[string]ontoCounts, cut, frac float64) error {
	if len(goTerms) == 0 || len(data.geneIDs) == 0 {
		return nil
	}
	root = strip(root, "<obo:", ">")

	sort.Strings(goTerms)
	m := mat.NewDense(len(data.geneIDs), len(goTerms), nil) // Assume all samples have same genes.
	for sample, name := range data.names {
		for col, term := range goTerms {
			counts, ok := ontoData[term]
			if !ok {
				continue
			}
			for row, geneID := range data.geneIDs {
				if counts.vector[sample].Bit(row) == 0 {
					continue
				}
				m.Set(row, col, data.counts[geneID][sample])
			}
		}

		path := fmt.Sprintf("%s_%s_%03d", name, root, depth)
		err := optimalTruncation(path, m, cut, frac)
		if err != nil {
			log.Println(err)
		}
		err = writeMatrix(path, data.geneIDs, goTerms, m)
		if err != nil {
			return err
		}

		m.Zero()
	}

	return nil
}

// https://arxiv.org/abs/1305.5870
func optimalTruncation(path string, m *mat.Dense, cut, frac float64) error {
	var svd mat.SVD
	ok := svd.Factorize(m, mat.SVDThin)
	if !ok {
		return fmt.Errorf("could not factorise %q", path)
	}
	sigma := svd.Values(nil)

	sum := make([]float64, len(sigma))
	floats.CumSum(sum, sigma)
	floats.Scale(1/sum[len(sum)-1], sum)
	rFrac := idxAbove(frac, sum)

	sigma = sigma[:idxBelow(cut, sigma)]
	rows, cols := m.Dims()
	t := tau(rows, cols, sigma)
	rOpt := idxBelow(t, sigma)

	if false {
		fmt.Printf("%s: %d[:%d]\n\t%v\n\t%v\n", path, len(sigma), rOpt, sigma, sigma[:rOpt])
	}

	p := plot.New()
	p.Title.Text = fmt.Sprintf("Singular Values\n%s", path)
	p.Y.Scale = logScale{}
	p.Y.Tick.Marker = logTicks{}
	sigXYs := sliceToXYs(sigma)
	if len(sigXYs) != 0 {
		values, err := plotter.NewLine(sigXYs)
		if err != nil {
			return err
		}
		threshOpt, err := plotter.NewLine(plotter.XYs{{X: 0, Y: t}, {X: sigXYs[len(sigXYs)-1].X, Y: t}})
		if err != nil {
			return err
		}
		threshOpt.Color = color.RGBA{B: 255, A: 255}
		p.Add(values, threshOpt)
		if rFrac < len(sigma) {
			threshFrac, err := plotter.NewLine(plotter.XYs{{X: 0, Y: sigma[rFrac]}, {X: sigXYs[len(sigXYs)-1].X, Y: sigma[rFrac]}})
			if err != nil {
				return err
			}
			threshFrac.Color = color.RGBA{R: 255, A: 255}
			p.Add(threshFrac)
		}
	}
	return p.Save(18*vg.Centimeter, 15*vg.Centimeter, filepath.Join("plots", path+".png"))
}

func idxAbove(thresh float64, s []float64) int {
	for i, v := range s {
		if v > thresh {
			return i
		}
	}
	return len(s)
}

func idxBelow(thresh float64, s []float64) int {
	for i, v := range s {
		if v < thresh {
			return i
		}
	}
	return len(s)
}

// https://arxiv.org/abs/1305.5870 Eq. 4.
func tau(rows, cols int, values []float64) float64 {
	if len(values) == 0 {
		return 0
	}
	reverseFloats(values)
	m := stat.Quantile(0.5, 1, values, nil)
	reverseFloats(values)
	return omega(rows, cols) * m
}

func reverseFloats(f []float64) {
	for i, j := 0, len(f)-1; i < j; i, j = i+1, j-1 {
		f[i], f[j] = f[j], f[i]
	}
}

// https://arxiv.org/abs/1305.5870 Eq. 5.
func omega(rows, cols int) float64 {
	beta := float64(rows) / float64(cols)
	beta2 := beta * beta
	return 0.56*beta2*beta - 0.95*beta2 + 1.82*beta + 1.43
}

func sliceToXYs(s []float64) plotter.XYs {
	xy := make(plotter.XYs, len(s))
	for i, v := range s {
		if v == 0 {
			return xy[:i]
		}
		xy[i] = plotter.XY{X: float64(i), Y: v}
	}
	return xy
}

type logScale struct{}

func (logScale) Normalize(min, max, x float64) float64 {
	min = math.Max(min, 1e-16)
	max = math.Max(max, 1e-16)
	x = math.Max(x, 1e-16)
	logMin := math.Log(min)
	return (math.Log(x) - logMin) / (math.Log(max) - logMin)
}

type logTicks struct{ powers int }

func (t logTicks) Ticks(min, max float64) []plot.Tick {
	min = math.Max(min, 1e-16)
	max = math.Max(max, 1e-16)
	if t.powers < 1 {
		t.powers = 1
	}

	val := math.Pow10(int(math.Log10(min)))
	max = math.Pow10(int(math.Ceil(math.Log10(max))))
	var ticks []plot.Tick
	for val < max {
		for i := 1; i < 10; i++ {
			if i == 1 {
				ticks = append(ticks, plot.Tick{Value: val, Label: strconv.FormatFloat(val, 'e', 0, 64)})
			}
			if t.powers != 1 {
				break
			}
			ticks = append(ticks, plot.Tick{Value: val * float64(i)})
		}
		val *= math.Pow10(t.powers)
	}
	ticks = append(ticks, plot.Tick{Value: val, Label: strconv.FormatFloat(val, 'e', 0, 64)})

	return ticks
}

func writeMatrix(path string, rows, cols []string, data *mat.Dense) (err error) {
	f, err := os.Create(filepath.Join("matrices", path+".tsv"))
	if err != nil {
		return err
	}
	defer func() {
		err = f.Close()
	}()

	_, err = f.Write([]byte{'\t'})
	if err != nil {
		return err
	}
	_, err = f.WriteString(strings.Join(stripSlice(cols, "<obo:", ">"), "\t"))
	if err != nil {
		return err
	}
	_, err = f.Write([]byte{'\n'})
	if err != nil {
		return err
	}

	for r, id := range rows {
		_, err = f.WriteString(id)
		if err != nil {
			return err
		}
		for c := range cols {
			_, err = fmt.Fprintf(f, "\t%v", data.At(r, c))
			if err != nil {
				return err
			}
		}
		_, err = f.Write([]byte{'\n'})
		if err != nil {
			return err
		}
	}
	return nil
}

func stripSlice(s []string, prefix, suffix string) []string {
	n := make([]string, len(s))
	for i, e := range s {
		n[i] = strip(e, prefix, suffix)
	}
	return n
}

func strip(s, prefix, suffix string) string {
	return strings.TrimSuffix(strings.TrimPrefix(s, prefix), suffix)
}

// countData holds count data for a set of named samples each with a collection
// of features.
type countData struct {
	// names is the names of the samples.
	names []string

	// counts holds a set of counts for
	// features keyed by the map key.
	// The length of each []float64
	// must match the length of names
	// and indexing into the []float64
	// reflects indexing into names.
	counts map[string][]float64

	// geneIDs and geneIdx are mappings
	// between internal gene index and
	// external gene identifier.
	geneIDs []string
	geneIdx map[string]int
}

// mappingCounts returns the count data held in the file at path.
func mappingCounts(path string) (*countData, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r, err := gzip.NewReader(f)
	if err != nil {
		return nil, err
	}

	c := csv.NewReader(r)
	c.Comma = '\t'
	c.Comment = '#'

	labels, err := c.Read()
	if err != nil {
		if err == io.EOF {
			return nil, io.ErrUnexpectedEOF
		}
		return nil, err
	}
	if labels[0] != "Geneid" {
		return nil, fmt.Errorf(`unexpected first column name: %q != "Geneid"`, labels[0])
	}
	samples := labels[1:]

	data := make(map[string][]float64, len(samples))
	geneIdx := make(map[string]int)
	var geneIDs []string

	c.ReuseRecord = true
	for {
		counts, err := c.Read()
		if err != nil {
			if err != io.EOF {
				return nil, err
			}
			break
		}
		geneid := counts[0]
		geneIdx[geneid] = len(geneIDs)
		geneIDs = append(geneIDs, geneid)
		for i, f := range counts[1:] {
			v, err := strconv.ParseFloat(f, 64)
			if err != nil {
				return nil, fmt.Errorf("error parsing value for %q in sample %q: %v", geneid, samples[i], err)
			}
			data[geneid] = append(data[geneid], v)
		}
	}

	return &countData{
		names:   samples,
		counts:  data,
		geneIDs: geneIDs,
		geneIdx: geneIdx,
	}, nil
}

// ontologyGraph returns the graph for the ontology stored in an OBO in OWL
// file. The namespaces are not expanded to full IRI namespaces.
func ontologyGraph(path string, lean bool) (*gogo.Graph, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r, err := gzip.NewReader(f)
	if err != nil {
		return nil, err
	}

	g := gogo.NewGraph()
	dec, err := owl.NewDecoder(r)
	if err != nil {
		return nil, err
	}
	for {
		s, err := dec.UnmarshalLocal()
		if err != nil {
			if err != io.EOF {
				return nil, err
			}
			break
		}

		if lean {
			// Filter on statements that are actually used. The filter
			// can be more restrictive since we only need "GO subclass GO"
			// and "GO hasOBONamespace literal" statements, but this is
			// good enough.
			switch s.Predicate.Value {
			// This list must include all predicates used in traversals
			// except <local:annotates> which comes from connectGeneIDsTo.
			case "<rdfs:subClassOf>", "<oboInOwl:hasOBONamespace>":
			default:
				continue
			}
		}

		g.AddStatement(s)
	}

	return g, nil
}

// connectGeneIDsTo adds the statements in path to the destination graph.
// The statements are expected to have local IRI namespaces and be in the
// following form:
//
//   <obo:GO_0000000> <local:annotates> <ensembl:ENSG00000000000> .
//
// Only ENSG identifiers that match the names in the counts map are added
// to the graph.
func connectGeneIDsTo(dst *gogo.Graph, path string, counts map[string][]float64) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	r, err := gzip.NewReader(f)
	if err != nil {
		return err
	}

	dec := rdf.NewDecoder(r)
	for {
		s, err := dec.Unmarshal()
		if err != nil {
			if err == io.EOF {
				return nil
			}
			return err
		}

		// Only keep annotations needed for the given counts.
		id := strip(s.Object.Value, "<ensembl:", ">")
		if _, ok := counts[id]; !ok {
			continue
		}

		s.Subject.UID = 0
		s.Predicate.UID = 0
		s.Object.UID = 0
		dst.AddStatement(s)
	}
}

type ontoCounts struct {
	vector []big.Int
}

// distributeCounts performs a breadth-first traversal from each of the leaf-most
// terms associated with each of the genes held by data, for each of the sample
// in data independently.
// Each gene ontology aspect is analysed separately since the aspects are not
// connected. The analyses are performed in parallel; the length of the
// returned slice will be the same as the number of roots passed in.
func distributeCounts(g *gogo.Graph, roots []rdf.Term, data *countData) []map[string]ontoCounts {
	ontoData := make([]map[string]ontoCounts, len(roots))
	for i := range ontoData {
		ontoData[i] = make(map[string]ontoCounts)
	}

	dfs := make([]traverse.DepthFirst, len(roots))
	for i := range dfs {
		dfs[i] = traverse.DepthFirst{Traverse: isSubClassOfGO}
	}
	for geneid, counts := range data.counts {
		var wg sync.WaitGroup
		for i, aspect := range leafiestFor(geneid, g, roots) {
			i := i
			aspect := aspect
			wg.Add(1)
			go func() {
				defer wg.Done()
				dfs[i].Reset()
				for _, l := range aspect {
					updateOntoData(ontoData[i], l, geneid, counts, data)
					dfs[i].Walk(g, l, func(n graph.Node) bool {
						t := n.(rdf.Term)
						updateOntoData(ontoData[i], t, geneid, counts, data)
						return false
					})
				}
			}()
		}
		wg.Wait()
	}

	return ontoData
}

func updateOntoData(ontoData map[string]ontoCounts, t rdf.Term, geneid string, counts []float64, data *countData) {
	dst, ok := ontoData[t.Value]
	if !ok {
		vector := make([]big.Int, len(counts))
		for _, v := range vector {
			v.SetBit(&v, len(data.geneIdx), 0)
		}
		dst.vector = vector
		ontoData[t.Value] = dst
	}
	for j, c := range counts {
		if c == 0 {
			continue
		}
		dst.vector[j].SetBit(&dst.vector[j], data.geneIdx[geneid], 1)
	}
}

// isSubClassOfGO is a traverse edge filter. It accepts statements where
//
//  any -- <rdfs:subClassOf> -> <obo:GO_*
//
// for out queries from a term.
func isSubClassOfGO(e graph.Edge) bool {
	return gogo.ConnectedByAny(e, func(s *rdf.Statement) bool {
		return s.Predicate.Value == "<rdfs:subClassOf>" &&
			strings.HasPrefix(s.Object.Value, "<obo:GO_")
	})
}

// leafiestFor return the leaf-most terms for gene from each of the ontology roots.
// The leaf sets are returned separated so that ontology count mutation can be
// performed concurrently without locking.
func leafiestFor(geneid string, g *gogo.Graph, roots []rdf.Term) [][]rdf.Term {
	leafiest := make([][]rdf.Term, len(roots))
	found := make([]bool, len(roots))
	var wg sync.WaitGroup
	for a, r := range roots {
		a := a
		r := r
		wg.Add(1)
		go func() {
			defer wg.Done()
			from, ok := g.TermFor("<ensembl:" + geneid + ">")
			if !ok {
				return
			}

			terms := g.Query(from).In(func(s *rdf.Statement) bool {
				return s.Predicate.Value == "<local:annotates>"
			}).Unique().Result()

			if len(terms) != 0 {
				found[a] = true
			}

			var depths []gogo.Descendant
			for _, q := range terms {
				ok, d := g.IsDescendantOf(r, q)
				if ok {
					depths = append(depths, gogo.Descendant{Term: q, Depth: d})
				}
			}
			sort.Sort(byDepth(depths))

			for i := 0; i < len(depths); i++ {
				a := depths[i]
				for j := i + 1; j < len(depths); {
					ok, _ := g.IsDescendantOf(a.Term, depths[j].Term)
					if ok {
						copy(depths[j:], depths[j+1:])
						depths = depths[:len(depths)-1]
					} else {
						j++
					}
				}
			}
			for _, d := range depths {
				leafiest[a] = append(leafiest[a], d.Term)
			}
		}()
	}
	wg.Wait()
	if !anyTrue(found) {
		log.Printf("no GO term found for %s", geneid)
	}
	return leafiest
}

func anyTrue(t []bool) bool {
	for _, v := range t {
		if v {
			return true
		}
	}
	return false
}

// byDepth sorts gogo.Descendents by depth, leafiest first.
type byDepth []gogo.Descendant

func (d byDepth) Len() int           { return len(d) }
func (d byDepth) Less(i, j int) bool { return d[i].Depth > d[j].Depth }
func (d byDepth) Swap(i, j int)      { d[i], d[j] = d[j], d[i] }

// nameSpaceOf returns the ontology aspect for the term t which is expected
// to be an <obo:GO_*> term.
func nameSpaceOf(t rdf.Term, in *gogo.Graph) string {
	ns := in.Query(t).Out(func(s *rdf.Statement) bool {
		return s.Predicate.Value == "<oboInOwl:hasOBONamespace>"
	}).Result()
	switch len(ns) {
	case 0:
		return "NA"
	case 1:
		text, _, kind, err := ns[0].Parts()
		if err != nil {
			panic(fmt.Errorf("invalid term in graph: %w", err))
		}
		if kind == rdf.Literal {
			return text
		}
		return t.Value
	default:
		return t.Value
	}
}

// walkDownSubClassesFrom performs a breadth-first enumeration of GO subclass
// terms in g starting from r, and calling fn for each term, including r.
func walkDownSubClassesFrom(r rdf.Term, g *gogo.Graph, fn func(root, term rdf.Term, depth int)) {
	bf := traverse.BreadthFirst{Traverse: goIsSubClassOf}
	bf.Walk(reverse{g}, r, func(n graph.Node, d int) bool {
		fn(r, n.(rdf.Term), d)
		return false
	})
}

// goIsSubClassOf is a traverse edge filter. It accepts statements where
//
//  <obo:GO_* <- <rdfs:subClassOf> -- any
//
// for in queries from a term.
func goIsSubClassOf(e graph.Edge) bool {
	return gogo.ConnectedByAny(e, func(s *rdf.Statement) bool {
		return s.Predicate.Value == "<rdfs:subClassOf>" &&
			strings.HasPrefix(s.Subject.Value, "<obo:GO_")
	})
}

// reverse implements the traverse.Graph reversing the direction of edges.
type reverse struct {
	*gogo.Graph
}

func (g reverse) From(id int64) graph.Nodes      { return g.Graph.To(id) }
func (g reverse) Edge(uid, vid int64) graph.Edge { return g.Graph.Edge(vid, uid) }
