package probability

import (
	"testing"
)

func TestSample(_ *testing.T) {
	Sample(NewUniform(0, 1), NewGenerator(0), 10)
}
