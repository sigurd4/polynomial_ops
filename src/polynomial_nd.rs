use core::{ops::{Mul, MulAssign, Add}, marker::Destruct};

use array__ops::{ArrayOps, ArrayNdOps, ArrayNd};

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
    ///             1 + y*2 + y*y*3 + x*(4 + y*5 + y*y*6) + x*x*(7 + y*8 + y*y*9)
    ///         );
    ///     }
    /// }
    /// ```
    fn evaluate_as_polynomial_nd(self, x: [X; N]) -> Y;
}

pub const fn sum_dims<const D: usize>(dims: [usize; D]) -> usize
{
    let mut s = 0;
    let mut d = 0;
    while d < D
    {
        s += dims[d];
        d += 1;
    }
    s
}

impl<C, X, A, const N: usize> /*const*/ PolynomialNd<X, <X as Mul<C>>::Output, N> for A
where
    A: ArrayNd<N, ItemNd = C> + /*~const*/ ArrayNdOps<N, C, {A::FLAT_LENGTH}>,
    C: /*~const*/ Destruct + /*~const*/ Into<<X as Mul<C>>::Output>,
    X: /*~const*/ Mul<C> + /*~const*/ MulAssign + /*~const*/ Mul<Output = X> + Copy,
    <X as Mul<C>>::Output: /*~const*/ Default + /*~const*/ Add<Output = <X as Mul<C>>::Output> + /*~const*/ Destruct,
    [(); sum_dims(A::DIMENSIONS)]:
{
    fn evaluate_as_polynomial_nd(self, x: [X; N]) -> <X as Mul<C>>::Output
    {
        let index_offset_in_xn: [usize; N] = {
            let mut n_accum = 0;
            A::DIMENSIONS.map2(const |d| {
                let i0 = n_accum;
                n_accum += d;
                i0
            })
        };
        let xn = {
            let mut xn = [None; sum_dims(A::DIMENSIONS)];

            let mut n = 0;
            while n < N
            {
                let mut x_accum = x[n];
                let mut i = 1;

                while i < A::DIMENSIONS[n]
                {
                    xn[index_offset_in_xn[n] + i] = Some(x_accum);
                    x_accum *= x[n];
                    i += 1;
                }
    
                n += 1;
            }

            xn
        };

        /*unsafe {core::intrinsics::const_eval_select((A::DIMENSIONS,), do_nothing, print)};
        unsafe {core::intrinsics::const_eval_select((A::DIMENSIONS.reduce(Add::add).unwrap_or_default(),), do_nothing, print)};
        unsafe {core::intrinsics::const_eval_select((&xn,), do_nothing, print)};*/

        self.enumerate_nd()
            .map_nd(|(indices, c)| {
                if let Some(x) = indices.zip(index_offset_in_xn)
                    .map2(|(i, i0)| xn[i0 + i])
                    .reduce(|a, b| match a
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
                }
            }).flatten_nd_array()
            .reduce(Add::add)
            .unwrap_or_default()
    }
}

/*const fn do_nothing<T>(_xn: T)
where
    T: ~const Destruct
{

}

fn print<T>(xn: T)
where
    T: Debug
{
    println!("{:?}", xn)
}*/