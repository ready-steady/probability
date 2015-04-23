package probability

import (
	"testing"

	"github.com/ready-steady/probability/uniform"
)

func TestSample(_ *testing.T) {
	uniform := uniform.New(0, 1)
	Sample(uniform, 10)
}
