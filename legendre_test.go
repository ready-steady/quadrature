package quadrature

import (
	"testing"

	"github.com/ready-steady/assert"
)

func TestLegendre(t *testing.T) {
	x, w := Legendre(42)
	assert.EqualWithin(x, fixtureLegendre42X, 1e-14, t)
	assert.EqualWithin(w, fixtureLegendre42W, 1e-14, t)
}
