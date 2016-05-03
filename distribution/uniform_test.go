package distribution

import (
	"testing"

	"github.com/ready-steady/assert"
)

func TestUniformCumulate(t *testing.T) {
	distribution := NewUniform(0.0, 1.0)
	x := []float64{-1.0, 0.5, 2.0}
	p := []float64{0.0, 0.5, 1.0}

	assert.Equal(Cumulate(distribution, x), p, t)
}

func TestUniformWeigh(t *testing.T) {
	distribution := NewUniform(0.0, 1.0)
	x := []float64{-1.0, 0.4269, 2.0}
	p := []float64{0.0, 1.0, 0.0}

	assert.Equal(Weigh(distribution, x), p, t)
}
