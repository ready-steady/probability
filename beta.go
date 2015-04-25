package probability

import (
	"github.com/ready-steady/special"
)

// Beta represents a beta distribution.
type Beta struct {
	α float64
	β float64
	a float64
	b float64
}

// NewBeta returns a beta distribution with α and β on [a, b].
func NewBeta(α, β, a, b float64) *Beta {
	return &Beta{α, β, a, b}
}

// CDF evaluates the CDF of the distribution.
func (b *Beta) CDF(x float64) float64 {
	return special.IncBeta((x-b.a)/(b.b-b.a), b.α, b.β, special.LogBeta(b.α, b.β))
}

// InvCDF evaluates the inverse CDF of the distribution.
func (b *Beta) InvCDF(p float64) float64 {
	return (b.b-b.a)*special.InvIncBeta(p, b.α, b.β, special.LogBeta(b.α, b.β)) + b.a
}
