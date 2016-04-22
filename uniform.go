package probability

// Uniform is a uniform distribution.
type Uniform struct {
	a float64
	b float64
}

// NewUniform returns a uniform distribution on [a, b].
func NewUniform(a, b float64) *Uniform {
	return &Uniform{a, b}
}

// Sample draws a sample.
func (self *Uniform) Sample(generator Generator) float64 {
	return (self.b-self.a)*generator.Float64() + self.a
}
