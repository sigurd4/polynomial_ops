use core::{ops::{Add, Mul, Sub, Neg, AddAssign}, alloc::Allocator, iter::Sum};

use num_traits::{Zero, One};
use num_identities_const::{OneConst, ZeroConst};
use array__ops::ArrayOps;

use super::*;

#[derive(Clone, Copy, Debug)]
pub struct ChebyshevPolynomial
{
    pub kind: usize,
    pub order: usize
}

impl ChebyshevPolynomial
{
    pub const fn new(kind: usize, order: usize) -> Self
    {
        Self {
            kind,
            order
        }
    }

    pub const fn new_of_first_kind(order: usize) -> Self
    {
        Self::new(1, order)
    }

    pub const fn new_of_second_kind(order: usize) -> Self
    {
        Self::new(2, order)
    }
}

#[cfg(feature = "std")]
impl<C, A> Into<Vec<C, A>> for ChebyshevPolynomial
where
    A: Allocator + Default + Clone,
    Vec<C, A>: IntoIterator<Item = C> + Polynomial<C> + FromIterator<C>,
    C: One + Zero + Neg<Output = C> + Mul<Output = C> + Sub<Output = C> + Sum + Clone
{
    fn into(self) -> Vec<C, A>
    {
        let mut t_prev: Vec<C, A> = [C::one()].into_iter().collect();

        if self.order == 0
        {
            return t_prev
        }

        let kind: C = (0..self.kind).map(|_| C::one()).sum();
        let mut t: Vec<C, A> = [C::zero(), kind].into_iter().collect();

        if self.order == 1
        {
            return t
        }

        let two = C::one() + C::one();
        for _ in 1..self.order
        {
            let mut t_prev_iter = t_prev.into_iter();
            let t_next = [-t_prev_iter.next().unwrap()].into_iter()
                .chain(t.clone().into_iter()
                    .zip(t_prev_iter.chain([C::zero(), C::zero()].into_iter()))
                    .map(|(tn, tn_prev)| two.clone() * tn - tn_prev)
                ).collect();
            t_prev = t;
            t = t_next;
        }

        t
    }
}

/*default impl<C, const N: usize> Into<[C; N]> for ChebyshevPolynomial
where
    [C; N]: Polynomial<<[C; N] as IntoIterator>::Item> + ArrayOps<C, N>,
    C: One + Zero + Copy
        + Add<Output = C> + Sub<Output = C> + Neg<Output = C> + AddAssign
        + Mul<Output = C>
{
    fn into(self) -> [C; N]
    {
        unsafe {
            chebyshev_array_with_one_and_zero_given_as(self.kind, self.order, C::zero(), C::one(), C::one() + C::one())
        }
    }
}*/

impl<C, const N: usize> const Into<Option<[C; N]>> for ChebyshevPolynomial
where
    [C; N]: Polynomial<<[C; N] as IntoIterator>::Item> + ArrayOps<C, N>,
    C: OneConst + ZeroConst + Copy
        + /*~const*/ Add<Output = C> + /*~const*/ Sub<Output = C> + /*~const*/ Neg<Output = C> + /*~const*/ AddAssign
        + /*~const*/ Mul<Output = C>
{
    fn into(self) -> Option<[C; N]>
    {
        unsafe {
            chebyshev_array_with_one_and_zero_given_as(self.kind, self.order, C::ZERO, C::ONE, C::ONE + C::ONE)
        }
    }
}

impl<C, const N: usize> const TryInto<[C; N]> for ChebyshevPolynomial
where
    [C; N]: Polynomial<<[C; N] as IntoIterator>::Item> + ArrayOps<C, N>,
    C: OneConst + ZeroConst + Copy
        + /*~const*/ Add<Output = C> + /*~const*/ Sub<Output = C> + /*~const*/ Neg<Output = C> + /*~const*/ AddAssign
        + /*~const*/ Mul<Output = C>
{
    type Error = Self;

    fn try_into(self) -> Result<[C; N], Self::Error>
    {
        Into::<Option<[C; N]>>::into(self).ok_or(self)
    }
}

/*const*/ unsafe fn chebyshev_array_with_one_and_zero_given_as<A, C, const N: usize>(kind: usize, order: usize, zero: C, one: C, two: C) -> Option<A>
where
    A: /*~const*/ ArrayOps<C, N> + Copy,
    C: Copy
        + /*~const*/ Add<Output = C> + /*~const*/ Sub<Output = C> + /*~const*/ Neg<Output = C> + /*~const*/ AddAssign
        + /*~const*/ Mul<Output = C>
{
    if order > N
    {
        return None
    }

    let mut t_prev: A = ArrayOps::fill(const |i| if i == 0 {one} else {zero});
    if order == 0
    {
        return Some(t_prev)
    }
    
    let mut kind_c = zero;
    let mut k = 0;
    while k < kind
    {
        kind_c += one;
        k += 1;
    }

    let mut t: A = ArrayOps::fill(const |i| if i == 1 {kind_c} else {zero});

    let mut k = 1;
    while k < order
    {
        let mut t_prev_iter = t_prev.into_iter();
        let mut t_iter = t.into_iter();
        let mut first = true;
        
        let t_next = ArrayOps::fill(/*const*/ |_| if first
            {
                first = false;
                -t_prev_iter.next().unwrap()
            }
            else
            {
                two * t_iter.next().unwrap() - t_prev_iter.next().unwrap()
            }
        );

        t_prev = t;
        t = t_next;
        k += 1;
    }

    Some(t)
}

impl<T> const Polynomial<T> for ChebyshevPolynomial
where
    T: OneConst + ZeroConst + Copy
        + /*~const*/ Add<T, Output = T> + /*~const*/ AddAssign
        + /*~const*/ Sub<T, Output = T>
        + /*~const*/ Mul<T, Output = T>
{
    type Y = T;

    fn evaluate_as_polynomial(self, x: T) -> Self::Y
    {
        unsafe {
            evaluate_chebyshev_given_one_and_zero(self, x, T::ZERO, T::ONE, T::ONE + T::ONE)
        }
    }
}

/*const*/ unsafe fn evaluate_chebyshev_given_one_and_zero<T>(p: ChebyshevPolynomial, x: T, zero: T, one: T, two: T) -> T
where
    T: Copy
        + /*~const*/ AddAssign
        + /*~const*/ Sub<T, Output = T>
        + /*~const*/ Mul<T, Output = T>
{
    let mut t_prev = one;
    if p.order == 0
    {
        return t_prev;
    }
    
    let mut kind_c = zero;
    let mut k = 0;
    while k < p.kind
    {
        kind_c += one;
        k += 1;
    }
    
    let mut t = x*kind_c;

    let mut k = 1;
    while k < p.order
    {
        let t_next = two*x*t - t_prev;
        t_prev = t;
        t = t_next;
        k += 1;
    }

    t
}