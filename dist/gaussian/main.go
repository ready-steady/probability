// Package gaussian provides algorithms for working with Gaussian distributions.
//
// https://en.wikipedia.org/wiki/Normal_distribution
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

// New returns a Gaussian distribution with mean μ and variance σ2.
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
	values := make([]float64, len(points))

	a, b := s.μ, math.Sqrt(2*s.σ2)

	for i, x := range points {
		values[i] = (1 + math.Erf((x-a)/b)) / 2
	}

	return values
}

// InvCDF evaluates the inverse CDF of the distribution.
func (s *Self) InvCDF(points []float64) []float64 {
	// Author: John Burkardt
	// Source: http://people.sc.fsu.edu/~jburkardt/c_src/asa241/asa241.html
	// Modified: December 27, 2004

	const (
		const1 = 0.180625
		const2 = 1.6
		split1 = 0.425
		split2 = 5
	)

	values := make([]float64, len(points))

	μ, σ := s.μ, math.Sqrt(s.σ2)
	inf := math.Inf(1)

	var q, r float64

	for i, p := range points {
		if p <= 0 {
			values[i] = -inf
			continue
		}
		if 1 <= p {
			values[i] = inf
			continue
		}

		q = p - 0.5

		if math.Abs(q) <= split1 {
			r = const1 - q*q
			values[i] = μ + σ*q*polynom(coefA, r)/polynom(coeffB, r)
			continue
		}

		if q < 0 {
			r = p
		} else {
			r = 1 - p
		}

		r = math.Sqrt(-math.Log(r))

		if r <= split2 {
			r -= const2
			values[i] = polynom(coefC, r) / polynom(coefD, r)
		} else {
			r -= split2
			values[i] = polynom(coefE, r) / polynom(coefF, r)
		}

		if q < 0 {
			values[i] = μ - σ*values[i]
		} else {
			values[i] = μ + σ*values[i]
		}
	}

	return values
}

func polynom(coef []float64, x float64) (value float64) {
	for i := 8 - 1; 0 <= i; i-- {
		value = value*x + coef[i]
	}

	return
}

var coefA = []float64{
	3.3871328727963666080,
	1.3314166789178437745e+2,
	1.9715909503065514427e+3,
	1.3731693765509461125e+4,
	4.5921953931549871457e+4,
	6.7265770927008700853e+4,
	3.3430575583588128105e+4,
	2.5090809287301226727e+3,
}

var coeffB = []float64{
	1.0,
	4.2313330701600911252e+1,
	6.8718700749205790830e+2,
	5.3941960214247511077e+3,
	2.1213794301586595867e+4,
	3.9307895800092710610e+4,
	2.8729085735721942674e+4,
	5.2264952788528545610e+3,
}

var coefC = []float64{
	1.42343711074968357734,
	4.63033784615654529590,
	5.76949722146069140550,
	3.64784832476320460504,
	1.27045825245236838258,
	2.41780725177450611770e-1,
	2.27238449892691845833e-2,
	7.74545014278341407640e-4,
}

var coefD = []float64{
	1.0,
	2.05319162663775882187,
	1.67638483018380384940,
	6.89767334985100004550e-1,
	1.48103976427480074590e-1,
	1.51986665636164571966e-2,
	5.47593808499534494600e-4,
	1.05075007164441684324e-9,
}

var coefE = []float64{
	6.65790464350110377720,
	5.46378491116411436990,
	1.78482653991729133580,
	2.96560571828504891230e-1,
	2.65321895265761230930e-2,
	1.24266094738807843860e-3,
	2.71155556874348757815e-5,
	2.01033439929228813265e-7,
}

var coefF = []float64{
	1.0,
	5.99832206555887937690e-1,
	1.36929880922735805310e-1,
	1.48753612908506148525e-2,
	7.86869131145613259100e-4,
	1.84631831751005468180e-5,
	1.42151175831644588870e-7,
	2.04426310338993978564e-15,
}
