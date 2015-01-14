// Package beta provides algorithms for working with beta distributions.
//
// https://en.wikipedia.org/wiki/Beta_distribution
package beta

import (
	"github.com/ready-steady/sfunc"
)

// Distribution represents a particular distribution from the family.
type Distribution struct {
	α float64
	β float64
	a float64
	b float64
}

// New returns a beta distribution with α and β on [a, b].
func New(α, β, a, b float64) *Distribution {
	return &Distribution{α, β, a, b}
}

// CDF evaluates the CDF of the distribution.
func (d *Distribution) CDF(x float64) float64 {
	return sfunc.IncBeta((x-d.a)/(d.b-d.a), d.α, d.β, sfunc.LogBeta(d.α, d.β))
}

// InvCDF evaluates the inverse CDF of the distribution.
func (d *Distribution) InvCDF(p float64) float64 {
	return (d.b-d.a)*sfunc.InvIncBeta(p, d.α, d.β, sfunc.LogBeta(d.α, d.β)) + d.a
}
