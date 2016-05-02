package distribution

import (
	"testing"
)

func TestSample(_ *testing.T) {
	Sample(NewUniform(0.0, 1.0), NewGenerator(0), 10)
}
