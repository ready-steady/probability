package distribution

import (
	"testing"

	"github.com/ready-steady/assert"
)

func TestUniformWeigh(t *testing.T) {
	distribution := NewUniform(0.0, 1.0)
	x := []float64{0.4269}
	p := []float64{1.0}

	assert.Equal(Weigh(distribution, x), p, t)
}
