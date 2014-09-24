// Package beta provides algorithms for working with beta distributions.
//
// https://en.wikipedia.org/wiki/Beta_distribution
package beta

import (
	"github.com/go-math/sfunc"
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
func (s *Self) CDF(x float64) float64 {
	return sfunc.IncBeta((x-s.a)/(s.b-s.a), s.α, s.β, sfunc.LogBeta(s.α, s.β))
}

// InvCDF evaluates the inverse CDF of the distribution.
func (s *Self) InvCDF(p float64) float64 {
	return (s.b-s.a)*sfunc.InvIncBeta(p, s.α, s.β, sfunc.LogBeta(s.α, s.β)) + s.a
}
