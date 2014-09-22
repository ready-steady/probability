// Package beta provides algorithms for working with beta distributions.
//
// https://en.wikipedia.org/wiki/Beta_distribution
package beta

import (
	"math"
)

// Self represents a particular distribution from the family.
type Self struct {
	α float64
	β float64
	a float64
	b float64
}

// New returns a beta distribution with α and β on [a, b].
func New(α, β, a, b float64) *Self {
	return &Self{α, β, a, b}
}

// CDF evaluates the CDF of the distribution.
func (s *Self) CDF(points []float64) []float64 {
	values := make([]float64, len(points))

	α, β, k, b := s.α, s.β, s.b-s.a, s.a

	for i, x := range points {
		values[i] = betaInc(α, β, (x-b)/k)
	}

	return values
}

func betaInc(α, β, x float64) float64 {
	// Author: John Burkardt
	// Source: http://people.sc.fsu.edu/~jburkardt/c_src/prob/prob.html

	const (
		max = 1000
		tol = 1.0e-07
	)

	if x <= 0 {
		return 0
	} else if 1 <= x {
		return 1
	}

	αx, βx := x, 1-x

	var flip bool
	if α < (α+β)*x {
		α, αx, β, βx = β, βx, α, αx
		flip = true
	}

	// https://en.wikipedia.org/wiki/Integration_by_reduction_formulae
	rx := αx
	n := int(β + βx*(α+β))
	if n != 0 {
		rx /= βx
	}

	value := 1.0

	// FIXME: Might not converge. Panic?
	for i, temp, psq, term := 1, β-1, α+β, 1.0; i <= max; {
		term = term * temp * rx / (α + float64(i))

		value += term
		temp = math.Abs(term)

		if temp <= tol && temp <= tol*value {
			break
		}

		i++
		n--

		if 0 <= n {
			temp = β - float64(i)
			if n == 0 {
				rx = αx
			}
		} else {
			temp = psq
			psq += 1
		}
	}

	// http://dlmf.nist.gov/8.17#v
	value = value * math.Exp(α*math.Log(αx)+(β-1)*math.Log(βx)) / (beta(α, β) * α)

	if flip {
		return 1 - value
	} else {
		return value
	}
}

func beta(x, y float64) float64 {
	z, _ := math.Lgamma(x + y)
	x, _ = math.Lgamma(x)
	y, _ = math.Lgamma(y)

	return math.Exp(x + y - z)
}
