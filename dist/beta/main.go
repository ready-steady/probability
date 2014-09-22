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

	var flip bool
	var i, ns int
	var cx, pp, psq, qq, rx, temp, term, value, xx float64

	psq = α + β

	if α < (α+β)*x {
		xx = 1 - x
		cx = x
		pp, qq = β, α
		flip = true
	} else {
		xx = x
		cx = 1 - x
		pp, qq = α, β
	}

	term = 1
	i = 1
	value = 1

	ns = int(qq + cx*(α+β))
	rx = xx / cx

	temp = qq - float64(i)
	if ns == 0 {
		rx = xx
	}

	// FIXME: Might not converge. Panic?
	for k := 0; k < max; k++ {
		term = term * temp * rx / (pp + float64(i))
		value = value + term
		temp = math.Abs(term)

		if temp <= tol && temp <= tol*value {
			break
		}

		i = i + 1
		ns = ns - 1

		if 0 <= ns {
			temp = qq - float64(i)
			if ns == 0 {
				rx = xx
			}
		} else {
			temp = psq
			psq = psq + 1
		}
	}

	value = value * math.Exp(pp*math.Log(xx)+(qq-1)*math.Log(cx)) / (beta(α, β) * pp)

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
