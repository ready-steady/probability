package distribution

import (
	"math/rand"
)

// Generator is a generator of random numbers.
type Generator interface {
	Float64() float64
	NormFloat64() float64
}

// NewGenerator returns a new generator.
func NewGenerator(seed int64) Generator {
	return rand.New(rand.NewSource(seed))
}
