// Package beta provides algorithms for working with beta distributions.
//
// https://en.wikipedia.org/wiki/Beta_distribution
package beta

import (
	"github.com/ready-steady/sfunc"
)

// Beta represents a particular distribution from the family.
type Beta struct {
	α float64
	β float64
	a float64
	b float64
}

// New returns a beta distribution with α and β on [a, b].
func New(α, β, a, b float64) *Beta {
	return &Beta{α, β, a, b}
}

// CDF evaluates the CDF of the distribution.
func (b *Beta) CDF(x float64) float64 {
	return sfunc.IncBeta((x-b.a)/(b.b-b.a), b.α, b.β, sfunc.LogBeta(b.α, b.β))
}

// InvCDF evaluates the inverse CDF of the distribution.
func (b *Beta) InvCDF(p float64) float64 {
	return (b.b-b.a)*sfunc.InvIncBeta(p, b.α, b.β, sfunc.LogBeta(b.α, b.β)) + b.a
}
