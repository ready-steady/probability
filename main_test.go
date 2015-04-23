package probability

import (
	"testing"

	"github.com/ready-steady/probability/generator"
	"github.com/ready-steady/probability/uniform"
)

func TestSample(_ *testing.T) {
	Sample(uniform.New(0, 1), generator.New(0), 10)
}
