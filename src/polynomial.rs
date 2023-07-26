use core::{ops::{Mul, Add, MulAssign, AddAssign}, marker::Destruct, alloc::Allocator};

#[const_trait]
pub trait Polynomial<X>: Sized
{
    type Y;

    /// Evaluates a polynomial
    /// 
    /// # Example
    /// 
    /// ```rust
    /// use polynomial_ops::*;
    /// 
    /// let p = [1, 2, 3];
    /// for x in 0..256
    /// {
    ///     assert_eq!(
    ///         p.evaluate_as_polynomial(x),
    ///         1 + x*2 + x*x*3
    ///     );
    /// }
    /// ```
    fn evaluate_as_polynomial(self, x: X) -> Self::Y;
}

impl<C, X> const Polynomial<X> for &[C]
where
    C: ~const Into<<X as Mul<C>>::Output> + Copy,
    X: ~const Mul<C> + ~const MulAssign + Copy,
    <X as Mul<C>>::Output: ~const Default + ~const Add<Output = <X as Mul<C>>::Output> + ~const AddAssign + Default
{
    type Y = <X as Mul<C>>::Output;

    fn evaluate_as_polynomial(self, x: X) -> Self::Y
    {
        unsafe {
            core::intrinsics::const_eval_select((self, x,), slice_polynomial::evaluate_const, slice_polynomial::evaluate_at_rt)
        }
    }
}

#[cfg(feature = "std")]
impl<C, X, A> Polynomial<X> for Vec<C, A>
where
    A: Allocator,
    C: Into<<X as Mul<C>>::Output>,
    X: Mul<C> + MulAssign + Copy,
    <X as Mul<C>>::Output: Default + Add<Output = <X as Mul<C>>::Output> + Default
{
    type Y = <X as Mul<C>>::Output;

    fn evaluate_as_polynomial(self, x: X) -> Self::Y
    {
        array_polynomial::evaluate_at_rt(self, x)
    }
}

impl<C, X, const N: usize> const Polynomial<X> for [C; N]
where
    C: ~const Destruct + ~const Into<<X as Mul<C>>::Output>,
    X: ~const Mul<C> + ~const MulAssign + Copy,
    <X as Mul<C>>::Output: ~const Default + ~const Add<Output = <X as Mul<C>>::Output> + ~const Destruct
{
    type Y = <X as Mul<C>>::Output;

    fn evaluate_as_polynomial(self, x: X) -> Self::Y
    {
        unsafe {
            core::intrinsics::const_eval_select((self, x,), array_polynomial::evaluate_const, array_polynomial::evaluate_at_rt)
        }
    }
}

mod slice_polynomial
{
    use core::ops::{Mul, MulAssign, Add, AddAssign};

    use super::array_polynomial;

    #[cfg(test)]
    #[test]
    fn test_evaluation_equality()
    {
        use crate::ChebyshevPolynomial;

        type T = i128;
        const ORDER: usize = 3;
        const N: usize = ORDER + 1;

        const CHEB: Result<[T; N], ChebyshevPolynomial> = ChebyshevPolynomial::new_of_first_kind(ORDER).try_into();

        for x in i8::MIN as i128..i8::MAX as i128
        {
            assert_eq!(
                evaluate_at_rt(&CHEB.unwrap(), x),
                evaluate_const(&CHEB.unwrap(), x)
            );
        }
    }

    pub const fn evaluate_const<C, X>(polynomial: &[C], x: X) -> <X as Mul<C>>::Output
    where
        C: ~const Into<<X as Mul<C>>::Output> + Copy,
        X: ~const Mul<C> + ~const MulAssign + Copy,
        <X as Mul<C>>::Output: ~const Default + ~const Add<Output = <X as Mul<C>>::Output> + ~const AddAssign
    {
        let mut y = match polynomial.get(0)
        {
            Some(c0) => (*c0).into(),
            None => return Default::default()
        };
        let mut xn = x;
        let ptr_range = polynomial[1..].as_ptr_range();
        let mut ptr = ptr_range.start;
        while unsafe {ptr.offset_from(ptr_range.end)} < 0
        {
            y += xn*unsafe {*ptr};
            xn *= x;
            ptr = unsafe {ptr.add(1)};
        }
        y
    }

    pub fn evaluate_at_rt<C, X>(polynomial: &[C], x: X) -> <X as Mul<C>>::Output
    where
        C: Into<<X as Mul<C>>::Output> + Copy,
        X: Mul<C> + MulAssign + Copy,
        <X as Mul<C>>::Output: Default + Add<Output = <X as Mul<C>>::Output>
    {
        array_polynomial::evaluate_at_rt(polynomial.into_iter().map(|&c| c), x)
    }
}

mod array_polynomial
{
    use core::{ops::{Mul, MulAssign, Add}, marker::Destruct};

    use array_trait::ArrayOps;

    #[cfg(test)]
    #[test]
    fn test_evaluation_equality()
    {
        use crate::ChebyshevPolynomial;

        type T = i128;
        const ORDER: usize = 3;
        const N: usize = ORDER + 1;

        const CHEB: Result<[T; N], ChebyshevPolynomial> = ChebyshevPolynomial::new_of_first_kind(ORDER).try_into();

        for x in i8::MIN as i128..i8::MAX as i128
        {
            assert_eq!(
                evaluate_at_rt(CHEB.unwrap(), x),
                evaluate_const(CHEB.unwrap(), x)
            );
        }
    }

    pub const fn evaluate_const<C, A, X, const N: usize>(polynomial: A, x: X) -> <X as Mul<C>>::Output
    where
        A: ~const ArrayOps<C, N>,
        C: ~const Destruct + ~const Into<<X as Mul<C>>::Output>,
        X: ~const Mul<C> + ~const MulAssign + Copy,
        <X as Mul<C>>::Output: ~const Default + ~const Add<Output = <X as Mul<C>>::Output> + ~const Destruct
    {
        let mut one = true;
        let mut xn = x;
        polynomial
            .map2(const |c| if one {
                one = false;
                c.into()
            }
            else
            {
                let xc = xn*c;
                xn *= x;
                xc
            })
            .reduce(Add::add)
            .unwrap_or(Default::default())
    }
    
    pub fn evaluate_at_rt<C, A, X>(polynomial: A, x: X) -> <X as Mul<C>>::Output
    where
        A: IntoIterator<Item = C>,
        C: Into<<X as Mul<C>>::Output>,
        X: Mul<C> + MulAssign + Copy,
        <X as Mul<C>>::Output: Default + Add<Output = <X as Mul<C>>::Output>
    {
        let mut one = true;
        let mut xn = x;
        polynomial
            .into_iter()
            .map(|c| if one {
                one = false;
                c.into()
            }
            else
            {
                let xc = xn*c;
                xn *= x;
                xc
            })
            .reduce(|a, b| a + b)
            .unwrap_or(Default::default())
    }
}