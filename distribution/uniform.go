package distribution

// Uniform is a uniform distribution.
type Uniform struct {
	a float64
	b float64
}

// NewUniform returns a uniform distribution on [a, b].
func NewUniform(a, b float64) *Uniform {
	return &Uniform{a, b}
}

// Cumulate evaluates the CDF.
func (self *Uniform) Cumulate(x float64) float64 {
	return (x - self.a) / (self.b - self.a)
}

// Invert evaluates the inverse of the CDF.
func (self *Uniform) Invert(p float64) float64 {
	return (self.b-self.a)*p + self.a
}

// Sample draws a sample.
func (self *Uniform) Sample(generator Generator) float64 {
	return self.Invert(generator.Float64())
}

// Weigh evaluates the PDF.
func (self *Uniform) Weigh(_ float64) float64 {
	return 1.0 / (self.b - self.a)
}
