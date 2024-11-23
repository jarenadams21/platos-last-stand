use std::f64::consts::PI;

/// Define a complex number struct
#[derive(Clone, Copy, Debug)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    /// Create a new complex number
    fn new(re: f64, im: f64) -> Self {
        Complex { re, im }
    }

    /// Compute the complex conjugate
    fn conj(&self) -> Self {
        Complex {
            re: self.re,
            im: -self.im,
        }
    }

    /// Compute the modulus (magnitude) of the complex number
    fn modulus(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }

    /// Multiply two complex numbers
    fn mul(self, other: Self) -> Self {
        Complex {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }

    /// Add two complex numbers
    fn add(self, other: Self) -> Self {
        Complex {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }

    /// Subtract two complex numbers
    fn sub(self, other: Self) -> Self {
        Complex {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }

    /// Multiply complex number by scalar
    fn mul_scalar(&self, scalar: f64) -> Self {
        Complex {
            re: self.re * scalar,
            im: self.im * scalar,
        }
    }

    /// Divide complex number by scalar
    fn div_scalar(&self, scalar: f64) -> Self {
        Complex {
            re: self.re / scalar,
            im: self.im / scalar,
        }
    }

    /// Create a complex number from polar coordinates
    fn from_polar(r: f64, theta: f64) -> Self {
        Complex {
            re: r * theta.cos(),
            im: r * theta.sin(),
        }
    }
}

impl std::ops::Neg for Complex {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Complex {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl std::ops::Add for Complex {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        self.add(other)
    }
}

impl std::ops::Sub for Complex {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self.sub(other)
    }
}

impl std::ops::Mul for Complex {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self.mul(other)
    }
}

impl std::ops::Mul<f64> for Complex {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self {
        self.mul_scalar(scalar)
    }
}

impl std::ops::Div<f64> for Complex {
    type Output = Self;
    fn div(self, scalar: f64) -> Self {
        self.div_scalar(scalar)
    }
}

/// A struct representing a vector of complex numbers
#[derive(Debug, Clone)]
struct Vector {
    data: Vec<Complex>,
}

impl Vector {
    /// Create a new vector of given size with all elements initialized to zero
    fn zeros(size: usize) -> Self {
        Vector {
            data: vec![Complex::new(0.0, 0.0); size],
        }
    }

    /// Create a vector from a slice of complex numbers
    fn from_slice(data: &[Complex]) -> Self {
        Vector {
            data: data.to_vec(),
        }
    }

    /// Get the dimension of the vector
    fn dim(&self) -> usize {
        self.data.len()
    }

    /// Indexing into the vector
    fn get(&self, index: usize) -> Complex {
        self.data[index]
    }

    /// Set value at index
    fn set(&mut self, index: usize, value: Complex) {
        self.data[index] = value;
    }

    /// Add another vector to this vector
    fn add(&self, other: &Vector) -> Self {
        let mut result = Vector::zeros(self.dim());
        for i in 0..self.dim() {
            result.data[i] = self.data[i] + other.data[i];
        }
        result
    }

    /// Multiply vector by a scalar
    fn mul_scalar(&self, scalar: Complex) -> Self {
        let mut result = Vector::zeros(self.dim());
        for i in 0..self.dim() {
            result.data[i] = self.data[i] * scalar;
        }
        result
    }

    /// Inner product with another vector (complex conjugate of self)
    fn inner_product(&self, other: &Vector) -> Complex {
        let mut sum = Complex::new(0.0, 0.0);
        for i in 0..self.dim() {
            sum = sum + self.get(i).conj() * other.get(i);
        }
        sum
    }

    /// Normalize the vector
    fn normalize(&self) -> Self {
        let norm = self.inner_product(self).re.sqrt();
        self.mul_scalar(Complex::new(1.0 / norm, 0.0))
    }
}

/// A struct representing a matrix of complex numbers
#[derive(Debug, Clone)]
struct Matrix {
    data: Vec<Vec<Complex>>,
    rows: usize,
    cols: usize,
}

impl Matrix {
    /// Create a new matrix with given dimensions, initialized to zero
    fn zeros(rows: usize, cols: usize) -> Self {
        Matrix {
            data: vec![vec![Complex::new(0.0, 0.0); cols]; rows],
            rows,
            cols,
        }
    }

    /// Create a matrix from a 2D array
    fn from_array(array: &[&[Complex]]) -> Self {
        let rows = array.len();
        let cols = if rows > 0 { array[0].len() } else { 0 };
        let mut data = Vec::with_capacity(rows);
        for row in array {
            data.push(row.to_vec());
        }
        Matrix { data, rows, cols }
    }

    /// Get element at (row, col)
    fn get(&self, row: usize, col: usize) -> Complex {
        self.data[row][col]
    }

    /// Set element at (row, col)
    fn set(&mut self, row: usize, col: usize, value: Complex) {
        self.data[row][col] = value;
    }

    /// Multiply two matrices
    fn mul(&self, other: &Matrix) -> Self {
        assert_eq!(self.cols, other.rows);
        let mut result = Matrix::zeros(self.rows, other.cols);
        for i in 0..self.rows {
            for j in 0..other.cols {
                let mut sum = Complex::new(0.0, 0.0);
                for k in 0..self.cols {
                    sum = sum + self.get(i, k) * other.get(k, j);
                }
                result.set(i, j, sum);
            }
        }
        result
    }

    /// Multiply matrix by a vector
    fn mul_vector(&self, vector: &Vector) -> Vector {
        assert_eq!(self.cols, vector.dim());
        let mut result = Vector::zeros(self.rows);
        for i in 0..self.rows {
            let mut sum = Complex::new(0.0, 0.0);
            for j in 0..self.cols {
                sum = sum + self.get(i, j) * vector.get(j);
            }
            result.set(i, sum);
        }
        result
    }

    /// Transpose the matrix
    fn transpose(&self) -> Self {
        let mut result = Matrix::zeros(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                result.set(j, i, self.get(i, j));
            }
        }
        result
    }

    /// Conjugate transpose of the matrix
    fn conj_transpose(&self) -> Self {
        let mut result = Matrix::zeros(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                result.set(j, i, self.get(i, j).conj());
            }
        }
        result
    }

    /// Compute the commutator of two matrices
    fn commutator(&self, other: &Matrix) -> Self {
        self.mul(other).sub(&other.mul(self))
    }

    /// Subtract another matrix from this matrix
    fn sub(&self, other: &Matrix) -> Self {
        assert_eq!(self.rows, other.rows);
        assert_eq!(self.cols, other.cols);
        let mut result = Matrix::zeros(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                result.set(i, j, self.get(i, j) - other.get(i, j));
            }
        }
        result
    }

    /// Add another matrix to this matrix
    fn add(&self, other: &Matrix) -> Self {
        assert_eq!(self.rows, other.rows);
        assert_eq!(self.cols, other.cols);
        let mut result = Matrix::zeros(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                result.set(i, j, self.get(i, j) + other.get(i, j));
            }
        }
        result
    }

    /// Multiply matrix by scalar
    fn mul_scalar(&self, scalar: Complex) -> Self {
        let mut result = Matrix::zeros(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                result.set(i, j, self.get(i, j) * scalar);
            }
        }
        result
    }
}

/// A struct representing a Hilbert space
struct HilbertSpace {
    dimension: usize,
    basis: Vec<Vector>,
}

impl HilbertSpace {
    /// Creates a new Hilbert space with the given dimension
    fn new(dimension: usize) -> Self {
        let mut basis = Vec::new();
        for i in 0..dimension {
            let mut vec = Vector::zeros(dimension);
            vec.set(i, Complex::new(1.0, 0.0));
            basis.push(vec);
        }
        HilbertSpace { dimension, basis }
    }

    /// Returns an orthonormal basis for the Hilbert space
    fn basis(&self) -> &Vec<Vector> {
        &self.basis
    }
}

/// Spin-1 operators
struct SpinOperators;

impl SpinOperators {
    /// Returns the spin-1 matrices Sx, Sy, Sz
    fn sx() -> Matrix {
        let factor = 1.0 / f64::sqrt(2.0);
        Matrix::from_array(&[
            &[Complex::new(0.0, 0.0), Complex::new(factor, 0.0), Complex::new(0.0, 0.0)],
            &[Complex::new(factor, 0.0), Complex::new(0.0, 0.0), Complex::new(factor, 0.0)],
            &[Complex::new(0.0, 0.0), Complex::new(factor, 0.0), Complex::new(0.0, 0.0)],
        ])
    }

    fn sy() -> Matrix {
        let factor = 1.0 / f64::sqrt(2.0);
        Matrix::from_array(&[
            &[Complex::new(0.0, 0.0), Complex::new(0.0, factor), Complex::new(0.0, 0.0)],
            &[Complex::new(0.0, -factor), Complex::new(0.0, 0.0), Complex::new(0.0, factor)],
            &[Complex::new(0.0, 0.0), Complex::new(0.0, -factor), Complex::new(0.0, 0.0)],
        ])
    }

    fn sz() -> Matrix {
        Matrix::from_array(&[
            &[Complex::new(1.0, 0.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0)],
            &[Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0)],
            &[Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)],
        ])
    }

    /// Compute the commutator [A, B]
    fn commutator(a: &Matrix, b: &Matrix) -> Matrix {
        a.mul(b).sub(&b.mul(a))
    }
}

/// Representing the quantum system with two interacting spin-1 particles
struct QuantumSystem {
    hilbert_space: HilbertSpace,
    state: Matrix, // Density matrix representation
}

impl QuantumSystem {
    /// Create a new quantum system with antisymmetric state
    fn new() -> Self {
        // Hilbert space of two spin-1 particles: dimension 3 x 3 = 9
        let hilbert_space = HilbertSpace::new(9);

        // Create antisymmetric (singlet-like) state for spin-1 particles
        let mut state = Matrix::zeros(9, 9);

        // Construct the antisymmetric state |Ψ^-⟩
        // For spin-1 particles, the antisymmetric subspace is 3-dimensional
        // For simplicity, we approximate an antisymmetric state

        // Basis vectors for single particle
        let basis_single = HilbertSpace::new(3).basis().clone();

        // Build the composite basis for two particles
        let mut composite_basis = Vec::new();
        for v_a in &basis_single {
            for v_b in &basis_single {
                let tensor_product = QuantumSystem::tensor_product(v_a, v_b);
                composite_basis.push(tensor_product);
            }
        }

        // Construct the antisymmetric state (simplified version)
        // |Ψ^-⟩ = (|1⟩| -1⟩ - | -1⟩|1⟩)/√2
        let idx_1 = 0; // Corresponds to |1⟩
        let idx_0 = 1; // Corresponds to |0⟩
        let idx_m1 = 2; // Corresponds to | -1⟩

        let mut psi_minus = Vector::zeros(9);
        let state_1_m1 = &composite_basis[idx_1 * 3 + idx_m1]; // |1⟩| -1⟩
        let state_m1_1 = &composite_basis[idx_m1 * 3 + idx_1]; // | -1⟩|1⟩

        for i in 0..9 {
            psi_minus.set(
                i,
                (state_1_m1.get(i) - state_m1_1.get(i)) * (1.0 / f64::sqrt(2.0)),
            );
        }

        // Create the density matrix ρ = |Ψ^-⟩⟨Ψ^-|
        for i in 0..9 {
            for j in 0..9 {
                let value = psi_minus.get(i).mul(psi_minus.get(j).conj());
                state.set(i, j, value);
            }
        }

        QuantumSystem {
            hilbert_space,
            state,
        }
    }

    /// Perform an operation on subsystem A and observe the effect on subsystem B
    fn apply_operator_on_a(&mut self, operator_a: &Matrix) {
        // Since we're dealing with two particles, we need to apply the operator on the composite system
        // The operator on the composite system is O_A ⊗ I_B
        let identity_b = Matrix::from_array(&[
            &[Complex::new(1.0, 0.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0)],
            &[Complex::new(0.0, 0.0), Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            &[Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
        ]);

        let operator_composite = QuantumSystem::tensor_product_matrix(operator_a, &identity_b);

        // Update the state: ρ' = O ρ O†
        let state_prime = operator_composite
            .mul(&self.state)
            .mul(&operator_composite.conj_transpose());

        self.state = state_prime;
    }

    /// Measure an observable on subsystem B
    fn measure_on_b(&self, operator_b: &Matrix) -> f64 {
        // The operator on the composite system is I_A ⊗ O_B
        let identity_a = Matrix::from_array(&[
            &[Complex::new(1.0, 0.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0)],
            &[Complex::new(0.0, 0.0), Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            &[Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
        ]);

        let operator_composite = QuantumSystem::tensor_product_matrix(&identity_a, operator_b);

        // Expectation value: ⟨O_B⟩ = Tr(ρ O_B)
        let mut trace = Complex::new(0.0, 0.0);
        for i in 0..self.state.rows {
            for j in 0..self.state.cols {
                trace = trace + self.state.get(i, j) * operator_composite.get(j, i);
            }
        }

        trace.re
    }

    /// Helper function to compute tensor product of two vectors
    fn tensor_product(v_a: &Vector, v_b: &Vector) -> Vector {
        let mut data = Vec::new();
        for a in &v_a.data {
            for b in &v_b.data {
                data.push(a.mul(*b));
            }
        }
        Vector { data }
    }

    /// Helper function to compute tensor product of two matrices
    fn tensor_product_matrix(a: &Matrix, b: &Matrix) -> Matrix {
        let rows = a.rows * b.rows;
        let cols = a.cols * b.cols;
        let mut data = vec![vec![Complex::new(0.0, 0.0); cols]; rows];

        for i in 0..a.rows {
            for j in 0..a.cols {
                for k in 0..b.rows {
                    for l in 0..b.cols {
                        let value = a.get(i, j) * b.get(k, l);
                        data[i * b.rows + k][j * b.cols + l] = value;
                    }
                }
            }
        }

        Matrix { data, rows, cols }
    }
}

fn main() {
    // Initialize the quantum system
    let mut system = QuantumSystem::new();

    // Define spin operators for subsystem A and B
    let sx = SpinOperators::sx();
    let sy = SpinOperators::sy();
    let sz = SpinOperators::sz();

    // Apply a rotation around the z-axis on subsystem A
    let theta = PI / 4.0; // Rotation angle

    // For spin-1 particles, the rotation operator around the z-axis is:
    // R_z(θ) = exp(-i θ S_z)
    // The eigenvalues of S_z are m = 1, 0, -1
    // So the rotation operator is a diagonal matrix with elements e^{-i θ m}

    let exp_i_theta = |m: f64| Complex::from_polar(1.0, -theta * m);

    let rotation_a = Matrix::from_array(&[
        &[exp_i_theta(1.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0)],
        &[Complex::new(0.0, 0.0), exp_i_theta(0.0), Complex::new(0.0, 0.0)],
        &[Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), exp_i_theta(-1.0)],
    ]);

    // Apply the rotation on subsystem A
    system.apply_operator_on_a(&rotation_a);

    // Measure Sx on subsystem B
    let sx_b = sx.clone();
    let expectation_sx_b = system.measure_on_b(&sx_b);

    println!(
        "Expectation value of Sx on subsystem B after rotation on A: {}",
        expectation_sx_b
    );
}
