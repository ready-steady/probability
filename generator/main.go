// Package generator provides a generator of random numbers.
package generator

import (
	"math/rand"
)

// Generator is a generator of random numbers.
type Generator interface {
	Float64() float64
	NormFloat64() float64
}

// New returns a new generator.
func New(seed int64) Generator {
	return rand.New(rand.NewSource(seed))
}
