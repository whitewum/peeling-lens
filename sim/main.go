// Command peeltraj runs degree-one peeling on a random k-uniform hypergraph
// near the peeling threshold and records the trajectory of the degree-one
// count Q_n together with the full degree histogram and two first-order
// approximations to Doob compensators.
//
// Each output row (one per recorded step) is
//
//	run n Q U2 S AQ AV U_0 ... U_40
//
// where n is the number of peeled edges, Q the number of degree-one
// vertices, U2 the number of degree-two vertices, S the number of remaining
// stubs, AQ the cumulative first-order drift of Q, AV the cumulative
// first-order drift of the weighted functional
//
//	V_n = sum_j f_j(a_c / a(n/M)) U_j(n),  f_j(b) = j b (1 - (1-b)^(j-1)),
//
// and U_j the number of vertices of degree j (degree 40 collects >= 40).
// AV is recorded only before the bottleneck time t*; after that it is NaN.
// Each run ends with a line "END run n d"; the run failed iff n < d.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"math"
	"math/rand"
	"os"
)

const maxDeg = 40

// critical returns x_c, q_c = 1 - e^{-x_c}, a_c = q_c^{k-1} and alpha_c.
func critical(k int) (x, q, ac, alpha float64) {
	lo, hi := 0.001, 10.0
	for i := 0; i < 100; i++ {
		m := (lo + hi) / 2
		if math.Expm1(m) > float64(k-1)*m {
			hi = m
		} else {
			lo = m
		}
	}
	x = (lo + hi) / 2
	q = -math.Expm1(-x)
	ac = math.Pow(q, float64(k-1))
	return x, q, ac, x / ac / float64(k)
}

// weights fills w[j] = f_j(b).
func weights(w []float64, b float64) {
	p := 1.0
	w[0], w[1] = 0, 0
	for j := 2; j <= maxDeg; j++ {
		p *= 1 - b
		w[j] = float64(j) * b * (1 - p)
	}
}

func main() {
	M := flag.Int("M", 100000, "number of vertices (cells)")
	k := flag.Int("k", 3, "edge size (cells per key)")
	r := flag.Float64("r", 0, "load offset: d/M = alpha_c + r/sqrt(M)")
	model := flag.String("model", "config", "config (iid endpoints) or plain (distinct cells per edge)")
	runs := flag.Int("runs", 20, "independent runs")
	stride := flag.Int("stride", 50, "record every stride peeling steps")
	seed := flag.Int64("seed", 1, "random seed")
	flag.Parse()
	if *model != "config" && *model != "plain" {
		fmt.Fprintln(os.Stderr, "model must be config or plain")
		os.Exit(2)
	}
	K := *k
	_, qc, ac, alpha := critical(K)
	d := int(float64(*M) * (alpha + *r/math.Sqrt(float64(*M))))
	mu := float64(K*d) / float64(*M)
	tstar := mu / float64(K) * (1 - math.Pow(qc, float64(K)))
	aOf := func(t float64) float64 {
		return math.Pow(math.Max((mu-float64(K)*t)/mu, 1e-12), float64(K-1)/float64(K))
	}

	w := bufio.NewWriter(os.Stdout)
	defer w.Flush()
	fmt.Fprintf(w, "# M=%d k=%d d=%d model=%s seed=%d\n", *M, K, d, *model, *seed)
	rng := rand.New(rand.NewSource(*seed))
	wNow := make([]float64, maxDeg+1)
	wNext := make([]float64, maxDeg+1)
	for run := 0; run < *runs; run++ {
		ends := make([]int32, d*K)
		deg := make([]int32, *M)
		for e := 0; e < d; e++ {
			for j := 0; j < K; j++ {
				v := int32(rng.Intn(*M))
				if *model == "plain" {
					for dup := true; dup; {
						dup = false
						for i := 0; i < j; i++ {
							if ends[e*K+i] == v {
								dup = true
								v = int32(rng.Intn(*M))
								break
							}
						}
					}
				}
				ends[e*K+j] = v
				deg[v]++
			}
		}
		// Incidence lists: vertex -> incident edges.
		start := make([]int32, *M+1)
		for _, v := range ends {
			start[v+1]++
		}
		for i := 0; i < *M; i++ {
			start[i+1] += start[i]
		}
		fill := append([]int32(nil), start[:*M]...)
		inc := make([]int32, d*K)
		for i, v := range ends {
			inc[fill[v]] = int32(i / K)
			fill[v]++
		}
		alive := make([]bool, d)
		for i := range alive {
			alive[i] = true
		}
		hist := make([]float64, maxDeg+1)
		Q, U2 := 0, 0
		stack := []int32{}
		for v := 0; v < *M; v++ {
			hist[min(int(deg[v]), maxDeg)]++
			switch deg[v] {
			case 1:
				Q++
				stack = append(stack, int32(v))
			case 2:
				U2++
			}
		}
		S := d * K
		AQ, AV := 0.0, 0.0
		n := 0
		emit := func() {
			av := AV
			if float64(n)/float64(*M) >= tstar {
				av = math.NaN()
			}
			fmt.Fprintf(w, "%d %d %d %d %d %.3f %.4f", run, n, Q, U2, S, AQ, av)
			for _, h := range hist {
				fmt.Fprintf(w, " %d", int(h))
			}
			fmt.Fprintln(w)
		}
		emit()
		for {
			var v int32 = -1
			for len(stack) > 0 {
				x := stack[len(stack)-1]
				stack = stack[:len(stack)-1]
				if deg[x] == 1 {
					v = x
					break
				}
			}
			if v < 0 {
				break
			}
			// One-step predictable drifts. The leaf stub is removed; the k-1
			// partner stubs are (to first order) uniform over the S-1 others.
			pool := float64(S - 1)
			if pool > 0 {
				AQ += -1 + float64(K-1)*float64(2*U2-(Q-1))/pool
				t := float64(n) / float64(*M)
				if t+1/float64(*M) < tstar {
					weights(wNow, ac/aOf(t))
					weights(wNext, ac/aOf(t+1/float64(*M)))
					// Histogram after the leaf moves from degree 1 to 0; f_0=f_1=0.
					var dw, hit float64
					for j := 2; j <= maxDeg; j++ {
						dw += (wNext[j] - wNow[j]) * hist[j]
						hit += float64(j) * hist[j] * (wNext[j-1] - wNext[j])
					}
					AV += dw + float64(K-1)/pool*hit
				}
			}
			var e int32 = -1
			for p := start[v]; p < start[v+1]; p++ {
				if alive[inc[p]] {
					e = inc[p]
					break
				}
			}
			alive[e] = false
			for j := 0; j < K; j++ {
				u := ends[int(e)*K+j]
				old := deg[u]
				deg[u]--
				nw := deg[u]
				hist[min(int(old), maxDeg)]--
				hist[min(int(nw), maxDeg)]++
				if old == 1 {
					Q--
				} else if old == 2 {
					U2--
				}
				if nw == 1 {
					Q++
					stack = append(stack, u)
				} else if nw == 2 {
					U2++
				}
			}
			S -= K
			n++
			if n%*stride == 0 {
				emit()
			}
		}
		if n%*stride != 0 {
			emit()
		}
		fmt.Fprintf(w, "END %d %d %d\n", run, n, d)
	}
}
