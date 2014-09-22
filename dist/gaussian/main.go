// Package gaussian provides algorithms for working with Gaussian distributions.
package gaussian

import (
	"math"
	"math/rand"
)

// Self represents a particular distribution from the family.
type Self struct {
	μ  float64
	σ2 float64
}

// New returns the distribution with mean μ and variance σ2.
func New(μ, σ2 float64) *Self {
	return &Self{μ, σ2}
}

// Sample draws samples from the distribution.
func (s *Self) Sample(count uint32) []float64 {
	points := make([]float64, count)

	μ, σ := s.μ, math.Sqrt(s.σ2)

	for i := range points {
		points[i] = μ + σ*rand.NormFloat64()
	}

	return points
}

// CDF evaluates the CDF of the distribution.
func (s *Self) CDF(points []float64) []float64 {
	// Author: John Burkardt
	// Source: http://people.sc.fsu.edu/~jburkardt/c_src/prob/prob.html

	const (
		a1  = 0.398942280444
		a2  = 0.399903438504
		a3  = 5.75885480458
		a4  = 29.8213557808
		a5  = 2.62433121679
		a6  = 48.6959930692
		a7  = 5.92885724438
		b0  = 0.398942280385
		b1  = 3.8052e-08
		b2  = 1.00000615302
		b3  = 3.98064794e-04
		b4  = 1.98615381364
		b5  = 0.151679116635
		b6  = 5.29330324926
		b7  = 4.8385912808
		b8  = 15.1508972451
		b9  = 0.742380924027
		b10 = 30.789933034
		b11 = 3.99019417011
	)

	cdf := make([]float64, len(points))

	μ, σ := s.μ, math.Sqrt(s.σ2)

	var absx, y, q float64

	for i, x := range points {
		x = (x - μ) / σ

		if x < 0 {
			absx = -x
		} else {
			absx = x
		}

		if absx <= 1.28 {
			y = 0.5 * x * x
			q = 0.5 - absx*(a1-a2*y/(y+a3-a4/(y+a5+a6/(y+a7))))
		} else if absx <= 12.7 {
			y = 0.5 * x * x
			q = math.Exp(-y) * b0 / (absx - b1 + b2/(absx+b3+b4/(absx-b5+b6/(absx+b7-b8/(absx+b9+b10/(absx+b11))))))
		} else {
			q = 0.0
		}

		if x < 0.0 {
			cdf[i] = q
		} else {
			cdf[i] = 1.0 - q
		}
	}

	return cdf
}
