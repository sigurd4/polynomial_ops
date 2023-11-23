use core::{mem::MaybeUninit, ops::{Mul, AddAssign}, alloc::Allocator};

use array_trait::ArrayOps;

use super::*;

#[const_trait]
pub trait ProductPolynomial
{
    type Output;

    /// Multiplies all polynomials together
    /// 
    /// # Example
    /// 
    /// ```rust
    /// #![feature(const_closures)]
    /// #![feature(const_trait_impl)]
    /// #![feature(const_mut_refs)]
    /// #![feature(generic_const_exprs)]
    /// 
    /// use array_trait::ArrayOps;
    /// use polynomial_ops::{Polynomial, ProductPolynomial};
    /// 
    /// // a = 1 + x
    /// const A: [u8; 2] = [1, 1];
    /// const X: [u8; 4] = [0, 1, 2, 3];
    /// const AX: [u8; 4] = X.map2(const |x| A.evaluate_as_polynomial(x));
    /// 
    /// // aa = a*a = (1 + x)^2 = 1 + 2x + x^2
    /// let aa = [A, A].product_polynomial();
    /// assert_eq!(aa, [1, 2, 1]);
    /// assert_eq!(X.map(|x| aa.evaluate_as_polynomial(x)), AX.map2(const |x| x*x));
    /// 
    /// // aaa = a*a*a = (1 + x)^3 = 1 + 3x + 3x^2 + x^3
    /// let aaa = [A, A, A].product_polynomial();
    /// assert_eq!(aaa, [1, 3, 3, 1]);
    /// assert_eq!(X.map(|x| aaa.evaluate_as_polynomial(x)), AX.map2(const |x| x*x*x));
    /// ```
    fn product_polynomial(self) -> Self::Output;
}

pub const fn polynomial_product_length(n: usize, m: usize) -> usize
{
    let k = n*m + 1 - m;
    if m > 1
    {
        k*m
    }
    else
    {
        k
    }
}

impl<T, const N: usize, const M: usize> /*const*/ ProductPolynomial for [[T; N]; M]
where
    T: /*~const*/ Default + /*~const*/ Mul<T, Output = T> + /*~const*/ AddAssign<T> + Copy,
    [(); polynomial_product_length(N, M)]:
{
    type Output = [T; polynomial_product_length(N, M)];

    fn product_polynomial(self) -> Self::Output
    {
        self.map2(/*const*/ |p| p.resize(/*const*/ |_| T::default()))
            .reduce(/*const*/ |a, b| {
                let len = polynomial_product_length(N, M);
                let mut y = [Default::default(); polynomial_product_length(N, M)];
                let mut k = 0;
                while k < len
                {
                    let k_next = k + 1;
        
                    let (mut i, mut j) = (k_next.saturating_sub(len), k_next.min(len));
                    let n = k_next.min(N);
        
                    while i < n
                    {
                        j -= 1;
        
                        y[k] += b[i]*a[j];
        
                        i += 1;
                    }
                    k = k_next;
                }
                y
            }).unwrap_or_else(const || unsafe {MaybeUninit::assume_init(MaybeUninit::uninit())})
    }
}

macro_rules! impl_product_polynomial_vec {
    ($(<{$($generics:tt)+}>)? for $type:ty $(where $($where:tt)+)?) => {
        impl<T, $($($generics)+)?> ProductPolynomial for $type
        where
            T: Default + Mul<T, Output = T> + AddAssign<T> + Copy,
            Vec<T>: MulPolynomial<Vec<T>, Output = Vec<T>>
            $(,$($where)+)?
        {
            type Output = Vec<T>;
    
            fn product_polynomial(self) -> Self::Output
            {
                self.into_iter()
                    .map(|p| p.to_vec())
                    .reduce(|a, b| a.mul_polynomial(b))
                    .unwrap_or_else(|| vec![])
            }
        }
    };
}

impl_product_polynomial_vec!(<{const N: usize, A}> for [Vec<T, A>; N] where A: Allocator);
impl_product_polynomial_vec!(<{const N: usize}> for [&[T]; N]);

impl_product_polynomial_vec!(<{const N: usize, A}> for &[Vec<T, A>; N] where A: Allocator);
impl_product_polynomial_vec!(<{const N: usize}> for &[&[T]; N]);

impl_product_polynomial_vec!(<{const N: usize}> for &[[T; N]]);
impl_product_polynomial_vec!(<{const N: usize}> for &[&[T; N]]);
impl_product_polynomial_vec!(for &[&[T]]);
impl_product_polynomial_vec!(<{A}> for &[Vec<T, A>] where A: Allocator);

impl_product_polynomial_vec!(<{const N: usize, A}> for Vec<[T; N], A> where A: Allocator);
impl_product_polynomial_vec!(<{const N: usize, A}> for Vec<&[T; N], A> where A: Allocator);
impl_product_polynomial_vec!(<{A}> for Vec<&[T], A> where A: Allocator);
impl_product_polynomial_vec!(<{A1, A2}> for Vec<Vec<T, A2>, A1> where A1: Allocator, A2: Allocator);

// M = 0 : 0
// M = 1 : N
// M = 2 : N + N - 1
// M = 3 : N + N - 1 + N - 1