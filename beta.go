package probability

import (
	"github.com/ready-steady/special"
)

// Beta is a beta distribution.
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

// Cumulate evaluates the CDF.
func (self *Beta) Cumulate(x float64) float64 {
	return special.IncBeta((x-self.a)/(self.b-self.a), self.α, self.β,
		special.LogBeta(self.α, self.β))
}

// Decumulate evaluates the inverse of the CDF.
func (self *Beta) Decumulate(p float64) float64 {
	return (self.b-self.a)*special.InvIncBeta(p, self.α, self.β,
		special.LogBeta(self.α, self.β)) + self.a
}
