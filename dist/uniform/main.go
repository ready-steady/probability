// Package uniform provides algorithms for working with uniform distributions.
package uniform

import (
	"math/rand"
)

type Self struct {
	a float64
	b float64
}

// New returns the (continuous) uniform distribution on [a, b].
func New(a, b float64) *Self {
	return &Self{a, b}
}

// Sample draws samples from the distribution.
func (s *Self) Sample(count uint32) []float64 {
	points := make([]float64, count)

	k, b := s.b-s.a, s.a

	for i := range points {
		// http://golang.org/src/pkg/math/rand/rand.go#L104
		points[i] = k*float64(rand.Int63n(1<<53))/(1<<53) + b
	}

	return points
}
