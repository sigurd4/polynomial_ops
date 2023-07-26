use core::{ops::{Mul, MulAssign, Add}, marker::Destruct};

use array_trait::{ArrayOps, ArrayNdOps, ArrayNd};

#[const_trait]
pub trait PolynomialNd<X, Y, const N: usize>: Sized
{
    /// Evaluates a multivariable polynomial
    /// 
    /// # Example
    /// 
    /// ```rust
    /// #![feature(generic_const_exprs)]
    /// 
    /// use polynomial_ops::*;
    /// 
    /// let p: [[u64; 3]; 3] = [
    ///     [1, 2, 3],
    ///     [4, 5, 6],
    ///     [7, 8, 9]
    /// ];
    /// for x in 0..256
    /// {
    ///     for y in 0..256
    ///     {
    ///         assert_eq!(
    ///             p.evaluate_as_polynomial_nd([x, y]),
    ///             1 + x*2 + x*x*3 + y*(4 + x*5 + x*x*6) + y*y*(7 + x*8 + x*x*9)
    ///         );
    ///     }
    /// }
    /// ```
    fn evaluate_as_polynomial_nd(self, x: [X; N]) -> Y;
}

impl<C, X, A, const N: usize> const PolynomialNd<X, <X as Mul<C>>::Output, N> for A
where
    A: ArrayNd<N, ItemNd = C> + ~const ArrayNdOps<N, C, {A::FLAT_LENGTH}>,
    C: ~const Destruct + ~const Into<<X as Mul<C>>::Output>,
    X: ~const Mul<C> + ~const MulAssign + ~const Mul<Output = X> + Copy,
    <X as Mul<C>>::Output: ~const Default + ~const Add<Output = <X as Mul<C>>::Output> + ~const Destruct,
    [(); A::DIMENSIONS.reduce(Add::add).unwrap_or_default()]:
{
    fn evaluate_as_polynomial_nd(self, x: [X; N]) -> <X as Mul<C>>::Output
    {
        let i0: [usize; N] = {
            let mut n_accum = 0;
            let mut dimensions_iter = A::DIMENSIONS.into_const_iter();
            ArrayOps::fill(const |_| {
                let i0 = n_accum;
                n_accum += dimensions_iter.next().unwrap();
                i0
            })
        };
        let xn = {
            let mut xn = [None; A::DIMENSIONS.reduce(Add::add).unwrap_or_default()];

            let mut n = 0;
            while n < N
            {
                if A::DIMENSIONS[n] < 2
                {
                    continue
                }
                let mut x_accum = x[n];
                let mut i = 1;

                while i < A::DIMENSIONS[n]
                {
                    xn[i0[n] + i] = Some(x_accum);
                    x_accum *= x[n];
                    i += 1;
                }
    
                n += 1;
            }

            xn
        };

        self.enumerate_nd()
            .map_nd(const |(i, c)| if let Some(x) = i.zip2(i0)
                .map2(const |(i, i0)| xn[i0 + i])
                .reduce(const |a, b| match a
                {
                    Some(a) => match b
                    {
                        Some(b) => Some(a*b),
                        None => Some(a)
                    },
                    None => b
                }).flatten()
            {
                x*c
            }
            else
            {
                c.into()
            }).flatten_nd_array()
            .reduce(Add::add)
            .unwrap_or_default()
    }
}