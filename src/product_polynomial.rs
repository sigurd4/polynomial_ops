use core::{mem::MaybeUninit, ops::{Mul, AddAssign}, alloc::Allocator};

use array_trait::ArrayOps;

use super::*;

#[const_trait]
pub trait ProductPolynomial
{
    type Output;

    fn product_polynomial(self) -> Self::Output;
}

impl<T, const N: usize, const M: usize> const ProductPolynomial for [[T; N]; M]
where
    T: ~const Default + ~const Mul<T, Output = T> + ~const AddAssign<T> + Copy,
    [(); (N*M + 1 - M)*M.min(1)]:
{
    type Output = [T; (N*M + 1 - M)*M.min(1)];

    fn product_polynomial(self) -> Self::Output
    {
        self.map2(const |p| p.resize(const |_| T::default()))
            .reduce(const |a, b| {
                let mut y = [const {Default::default()}; (N*M + 1 - M)*M.min(1)];
                let mut k = 0;
                while k < (N*M + 1 - M)*M.min(1)
                {
                    let k_next = k + 1;
        
                    let (mut i, mut j) = (k_next.saturating_sub((N*M + 1 - M)*M.min(1)), k_next.min((N*M + 1 - M)*M.min(1)));
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