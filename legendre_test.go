package quadrature

import (
	"testing"

	"github.com/ready-steady/assert"
)

func TestLegendre(t *testing.T) {
	x, w := Legendre(42, -1.0, 1.0)
	assert.Close(x, fixtureLegendre42X, 1e-14, t)
	assert.Close(w, fixtureLegendre42W, 1e-14, t)
}
