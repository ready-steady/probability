// Package sample provides algorithms for drawing samples from various
// probability distributions.
package sample

import (
	"math/rand"
)

// Seed sets the random number generator to a deterministic state.
func Seed(seed int64) {
	rand.Seed(seed)
}

// Uniform draws samples from the uniform distribution on [0, 1).
func Uniform(count uint32) []float64 {
	points := make([]float64, count)

	for i := range points {
		// http://golang.org/src/pkg/math/rand/rand.go#L104
		points[i] = float64(rand.Int63n(1<<53)) / (1 << 53)
	}

	return points
}
