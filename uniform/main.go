// Package uniform provides algorithms for working with uniform distributions.
//
// https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
package uniform

import (
	"math/rand"
)

// Distribution represents a particular distribution from the family.
type Distribution struct {
	a float64
	b float64
}

// New returns a uniform distribution on [a, b].
func New(a, b float64) *Distribution {
	return &Distribution{a, b}
}

// Sample draws a sample from the distribution.
func (d *Distribution) Sample() float64 {
	// http://golang.org/src/pkg/math/rand/rand.go#L104
	return (d.b-d.a)*float64(rand.Int63n(1<<53))/(1<<53) + d.a
}
