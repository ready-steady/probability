// Package uniform provides algorithms for working with uniform distributions.
//
// https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
package uniform

import (
	"github.com/ready-steady/probability/generator"
)

// Uniform represents a uniform distribution.
type Uniform struct {
	a float64
	b float64
}

// New returns a uniform distribution on [a, b].
func New(a, b float64) *Uniform {
	return &Uniform{a, b}
}

// Sample draws a sample from the distribution.
func (u *Uniform) Sample(generator generator.Generator) float64 {
	return (u.b-u.a)*generator.Float64() + u.a
}
