package distribution

import (
	"math"

	"github.com/ready-steady/special"
)

// Beta is a beta distribution.
type Beta struct {
	α float64
	β float64
	a float64
	b float64

	lnBeta float64
}

// NewBeta returns a beta distribution with α and β on [a, b].
func NewBeta(α, β, a, b float64) *Beta {
	return &Beta{α, β, a, b, special.LnBeta(α, β)}
}

// Cumulate evaluates the CDF.
func (self *Beta) Cumulate(x float64) float64 {
	return special.IncBeta((x-self.a)/(self.b-self.a), self.α, self.β, self.lnBeta)
}

// Invert evaluates the inverse of the CDF.
func (self *Beta) Invert(p float64) float64 {
	return (self.b-self.a)*special.InvIncBeta(p, self.α, self.β, self.lnBeta) + self.a
}

// Sample draws a sample.
func (self *Beta) Sample(generator Generator) float64 {
	return self.Invert(generator.Float64())
}

// Weigh evaluates the PDF.
func (self *Beta) Weigh(x float64) float64 {
	scale := self.b - self.a
	x = (x - self.a) / scale
	return math.Exp((self.α-1.0)*math.Log(x)+(self.β-1.0)*math.Log(1.0-x)-self.lnBeta) / scale
}
