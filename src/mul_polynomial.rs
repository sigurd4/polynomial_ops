use core::{alloc::Allocator, ops::{Mul, Add, AddAssign}};

#[const_trait]
pub trait MulPolynomial<Rhs>
{
    type Output;

    /// Multiplies two polynomials, and yields a new, larger polynomial
    /// 
    /// # Example
    /// 
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use polynomial_ops::*;
    /// 
    /// assert_eq!(
    ///     [1.0, 1.0].mul_polynomial([1.0, 1.0]),
    ///     [1.0, 2.0, 1.0]
    /// );
    /// assert_eq!(
    ///     [1.0, -1.0].mul_polynomial([1.0, -1.0]),
    ///     [1.0, -2.0, 1.0]
    /// );
    /// assert_eq!(
    ///     [1.0, 1.0].mul_polynomial([1.0, -1.0]),
    ///     [1.0, 0.0, -1.0]
    /// );
    /// ```
    fn mul_polynomial(self, rhs: Rhs) -> Self::Output;
}

impl<C1, C2, const N1: usize, const N2: usize> const MulPolynomial<[C2; N2]> for [C1; N1]
where
    C1: /*~const*/ Mul<C2> + Copy,
    C2: Copy,
    <C1 as Mul<C2>>::Output: /*~const*/ AddAssign<<C1 as Mul<C2>>::Output> + /*~const*/ Default + Copy,
    [(); N1 + N2 - 1]:
{
    type Output = [<C1 as Mul<C2>>::Output; N1 + N2 - 1];

    fn mul_polynomial(self, rhs: [C2; N2]) -> Self::Output
    {
        let mut y = [Default::default(); N1 + N2 - 1];
        let mut k = 0;
        while k < N1 + N2 - 1
        {
            let k_next = k + 1;

            let (mut i, mut j) = (k_next.saturating_sub(N2), k_next.min(N2));
            let n = k_next.min(N1);

            while i < n
            {
                j -= 1;

                y[k] += self[i]*rhs[j];

                i += 1;
            }
            k = k_next;
        }
        y
    }
}

#[cfg(feature = "std")]
impl<C1, C2> MulPolynomial<&[C2]> for &[C1]
where
    C1: Mul<C2> + Copy,
    C2: Copy,
    <C1 as Mul<C2>>::Output: Add<<C1 as Mul<C2>>::Output, Output = <C1 as Mul<C2>>::Output> + Default
{
    type Output = Vec<<C1 as Mul<C2>>::Output>;

    fn mul_polynomial(self, rhs: &[C2]) -> Self::Output
    {
        let self_len = self.len();
        let rhs_len = rhs.len();
        let len = self_len + rhs_len - 1;
        (1..=len)
            .map(|k| (k.saturating_sub(rhs_len)..k.min(self_len))
                .zip((k.saturating_sub(self_len)..k.min(rhs_len)).rev())
                .map(|(i, j)| self[i]*rhs[j])
                .reduce(|a, b| a + b)
                .unwrap_or_default()
            ).collect()
    }
}

#[cfg(feature = "std")]
impl<C1, A1, C2, A2> MulPolynomial<Vec<C2, A2>> for Vec<C1, A1>
where
    A1: Allocator,
    A2: Allocator,
    C1: Mul<C2> + Copy,
    C2: Copy,
    <C1 as Mul<C2>>::Output: Add<<C1 as Mul<C2>>::Output, Output = <C1 as Mul<C2>>::Output> + Default
{
    type Output = Vec<<C1 as Mul<C2>>::Output>;

    fn mul_polynomial(self, rhs: Vec<C2, A2>) -> Self::Output
    {
        self.as_slice().mul_polynomial(rhs.as_slice())
    }
}

#[cfg(feature = "std")]
impl<C1, A1, C2> MulPolynomial<&[C2]> for Vec<C1, A1>
where
    A1: Allocator,
    C1: Mul<C2> + Copy,
    C2: Copy,
    <C1 as Mul<C2>>::Output: Add<<C1 as Mul<C2>>::Output, Output = <C1 as Mul<C2>>::Output> + Default
{
    type Output = Vec<<C1 as Mul<C2>>::Output>;

    fn mul_polynomial(self, rhs: &[C2]) -> Self::Output
    {
        self.as_slice().mul_polynomial(rhs)
    }
}

#[cfg(feature = "std")]
impl<C1, C2, A2> MulPolynomial<Vec<C2, A2>> for &[C1]
where
    A2: Allocator,
    C1: Mul<C2> + Copy,
    C2: Copy,
    <C1 as Mul<C2>>::Output: Add<<C1 as Mul<C2>>::Output, Output = <C1 as Mul<C2>>::Output> + Default
{
    type Output = Vec<<C1 as Mul<C2>>::Output>;

    fn mul_polynomial(self, rhs: Vec<C2, A2>) -> Self::Output
    {
        self.mul_polynomial(rhs.as_slice())
    }
}



#[cfg(feature = "std")]
impl<C1, C2, const N: usize> MulPolynomial<&[C2]> for [C1; N]
where
    C1: Mul<C2> + Copy,
    C2: Copy,
    <C1 as Mul<C2>>::Output: Add<<C1 as Mul<C2>>::Output, Output = <C1 as Mul<C2>>::Output> + Default
{
    type Output = Vec<<C1 as Mul<C2>>::Output>;

    fn mul_polynomial(self, rhs: &[C2]) -> Self::Output
    {
        self.as_slice().mul_polynomial(rhs)
    }
}

#[cfg(feature = "std")]
impl<C1, C2, A, const N: usize> MulPolynomial<Vec<C2, A>> for [C1; N]
where
    A: Allocator,
    C1: Mul<C2> + Copy,
    C2: Copy,
    <C1 as Mul<C2>>::Output: Add<<C1 as Mul<C2>>::Output, Output = <C1 as Mul<C2>>::Output> + Default
{
    type Output = Vec<<C1 as Mul<C2>>::Output>;

    fn mul_polynomial(self, rhs: Vec<C2, A>) -> Self::Output
    {
        self.as_slice().mul_polynomial(rhs)
    }
}

#[cfg(feature = "std")]
impl<C1, C2, const N: usize> MulPolynomial<[C2; N]> for &[C1]
where
    C1: Mul<C2> + Copy,
    C2: Copy,
    <C1 as Mul<C2>>::Output: Add<<C1 as Mul<C2>>::Output, Output = <C1 as Mul<C2>>::Output> + Default
{
    type Output = Vec<<C1 as Mul<C2>>::Output>;

    fn mul_polynomial(self, rhs: [C2; N]) -> Self::Output
    {
        self.mul_polynomial(rhs.as_slice())
    }
}

#[cfg(feature = "std")]
impl<C1, C2, A, const N: usize> MulPolynomial<[C2; N]> for Vec<C1, A>
where
    A: Allocator,
    C1: Mul<C2> + Copy,
    C2: Copy,
    <C1 as Mul<C2>>::Output: Add<<C1 as Mul<C2>>::Output, Output = <C1 as Mul<C2>>::Output> + Default
{
    type Output = Vec<<C1 as Mul<C2>>::Output>;

    fn mul_polynomial(self, rhs: [C2; N]) -> Self::Output
    {
        self.mul_polynomial(rhs.as_slice())
    }
}