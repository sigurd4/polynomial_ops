use std::{ops::{Mul, AddAssign, Add, MulAssign, IndexMut, Index}, marker::Destruct, borrow::{BorrowMut, Borrow}};

use array_trait::{Array, ArrayOps};
use currying::{RCurry};
use fn_zip::FnZip;
use num_identities_const::{OneConst, ZeroConst};

use super::*;

#[derive(Clone, Copy, Debug)]
pub struct ArrayPolynomial<C, const N: usize>(pub [C; N]);

impl<C, const N: usize> ArrayPolynomial<C, N>
{
    pub const fn from(array: [C; N]) -> Self
    {
        Self(array)
    }
    pub const fn from_ref(array: &[C; N]) -> &Self
    {
        unsafe {core::mem::transmute(array)}
    }
    pub const fn from_mut(array: &mut [C; N]) -> &mut Self
    {
        unsafe {core::mem::transmute(array)}
    }

    pub const fn into_array(self) -> [C; N]
    {
        private::array_polynomial_into_inner(self)
    }
    pub const fn get_array(&self) -> &[C; N]
    {
        &self.0
    }
    pub const fn get_array_mut(&mut self) -> &mut [C; N]
    {
        &mut self.0
    }

    pub const fn evaluate<X>(self, x: X) -> <X as Mul<C>>::Output
    where
        C: ~const Destruct + ~const Into<<X as Mul<C>>::Output>,
        X: ~const Mul<C> + ~const MulAssign + Copy,
        <X as Mul<C>>::Output: ~const Default + ~const Add<Output = <X as Mul<C>>::Output> + ~const Destruct
    {
        unsafe {
            core::intrinsics::const_eval_select((self, x,), Self::evaluate_const, Self::evaluate_at_rt)
        }
    }

    pub const fn mul_each<M>(self, rhs: ArrayPolynomial<M, N>) -> ArrayPolynomial<<C as Mul<M>>::Output, N>
    where
        C: ~const Mul<M>,
        M: ~const Destruct
    {
        unsafe {
            core::intrinsics::const_eval_select((self, rhs,), Self::mul_each_const, Self::mul_each_at_rt)
        }
    }

    pub const fn mul_every<M>(self, rhs: M) -> ArrayPolynomial<<C as Mul<M>>::Output, N>
    where
        C: ~const Mul<M>,
        M: Copy
    {
        self.map2(&Mul::mul.rcurry(rhs))
    }


    const fn evaluate_const<X>(self, x: X) -> <X as Mul<C>>::Output
    where
        C: ~const Destruct + ~const Into<<X as Mul<C>>::Output>,
        X: ~const Mul<C> + ~const MulAssign + Copy,
        <X as Mul<C>>::Output: ~const Default + ~const Add<Output = <X as Mul<C>>::Output> + ~const Destruct
    {
        let mut one = true;
        let mut xn = x;
        self.0
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
    
    fn evaluate_at_rt<X>(self, x: X) -> <X as Mul<C>>::Output
    where
        C: Into<<X as Mul<C>>::Output>,
        X: Mul<C> + MulAssign + Copy,
        <X as Mul<C>>::Output: Default + Add<Output = <X as Mul<C>>::Output>
    {
        let mut one = true;
        let mut xn = x;
        self.0
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

    const fn mul_each_const<M>(self, rhs: ArrayPolynomial<M, N>) -> ArrayPolynomial<<C as Mul<M>>::Output, N>
    where
        C: ~const Mul<M>
    {
        let mut iter_self = self.into_const_iter();
        let mut iter_rhs = rhs.into_const_iter();
        ArrayPolynomial::fill(const |_| iter_self.next().unwrap()*iter_rhs.next().unwrap())
    }

    fn mul_each_at_rt<M>(self, rhs: ArrayPolynomial<M, N>) -> ArrayPolynomial<<C as Mul<M>>::Output, N>
    where
        C: Mul<M>
    {
        self.zip2(rhs).map2(|(c, m)| c*m)
    }

    const fn mul_polynomial<T, const M: usize>(self, rhs: ArrayPolynomial<T, M>)
        -> ArrayPolynomial<<C as Mul<T>>::Output, {N + M}>
    where
        C: ~const Mul<T> + Copy,
        T: Copy,
        <C as Mul<T>>::Output: ~const Default + ~const AddAssign,
        [(); N + M]:
    {
        let mut array = [const {Default::default()}; N + M];

        let mut i = 0;
        while i < N
        {
            let mut j = 0;
            while j < M
            {
                array[i + j] += self.0[i]*rhs.0[j];
                j += 1;
            }
            i += 1;
        }

        ArrayPolynomial(array)
    }
}
impl<C, T, const N: usize, const M: usize> const Mul<ArrayPolynomial<T, M>> for ArrayPolynomial<C, N>
where
    C: ~const Mul<T> + Copy,
    T: Copy,
    <C as Mul<T>>::Output: ~const Default + ~const AddAssign,
    [(); N + M]:
{
    type Output = ArrayPolynomial<<C as Mul<T>>::Output, {N + M}>;

    fn mul(self, rhs: ArrayPolynomial<T, M>) -> Self::Output
    {
        self.mul_polynomial(rhs)
    }
}

impl<C, const N: usize> const From<[C; N]> for ArrayPolynomial<C, N>
{
    fn from(value: [C; N]) -> Self
    {
        Self::from(value)
    }
}
impl<'a, C, const N: usize> const From<&'a [C; N]> for &'a ArrayPolynomial<C, N>
{
    fn from(value: &'a [C; N]) -> Self
    {
        ArrayPolynomial::from_ref(value)
    }
}
impl<'a, C, const N: usize> const From<&'a mut [C; N]> for &'a mut ArrayPolynomial<C, N>
{
    fn from(value: &'a mut [C; N]) -> Self
    {
        ArrayPolynomial::from_mut(value)
    }
}

impl<C, const N: usize> const Into<[C; N]> for ArrayPolynomial<C, N>
{
    fn into(self) -> [C; N]
    {
        self.into_array()
    }
}
impl<'a, C, const N: usize> const Into<&'a [C; N]> for &'a ArrayPolynomial<C, N>
{
    fn into(self) -> &'a [C; N]
    {
        self.get_array()
    }
}
impl<'a, C, const N: usize> const Into<&'a mut [C; N]> for &'a mut ArrayPolynomial<C, N>
{
    fn into(self) -> &'a mut [C; N]
    {
        self.get_array_mut()
    }
}

impl<C, const N: usize, I> const Index<I> for ArrayPolynomial<C, N>
where
    [C; N]: ~const Index<I>
{
    type Output = <[C; N] as Index<I>>::Output;

    fn index(&self, index: I) -> &Self::Output
    {
        self.get_array().index(index)
    }
}

impl<C, const N: usize, I> const IndexMut<I> for ArrayPolynomial<C, N>
where
    [C; N]: ~const IndexMut<I>
{
    fn index_mut(&mut self, index: I) -> &mut Self::Output
    {
        self.get_array_mut().index_mut(index)
    }
}

impl<C, const N: usize> const IntoIterator for ArrayPolynomial<C, N>
where
    [C; N]: ~const IntoIterator
{
    type Item = <[C; N] as IntoIterator>::Item;
    type IntoIter = <[C; N] as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter
    {
        self.into_array().into_iter()
    }
}

impl<C, const N: usize> const Borrow<[C]> for ArrayPolynomial<C, N>
where
    [C; N]: ~const Borrow<[C]>
{
    fn borrow(&self) -> &[C]
    {
        self.get_array().borrow()
    }
}
impl<C, const N: usize> const BorrowMut<[C]> for ArrayPolynomial<C, N>
where
    [C; N]: ~const BorrowMut<[C]>
{
    fn borrow_mut(&mut self) -> &mut [C]
    {
        self.get_array_mut().borrow_mut()
    }
}

impl<C, const N: usize> const Borrow<[C; N]> for ArrayPolynomial<C, N>
where
    [C; N]: ~const Borrow<[C; N]>
{
    fn borrow(&self) -> &[C; N]
    {
        self.get_array().borrow()
    }
}
impl<C, const N: usize> const BorrowMut<[C; N]> for ArrayPolynomial<C, N>
where
    [C; N]: ~const BorrowMut<[C; N]>
{
    fn borrow_mut(&mut self) -> &mut [C; N]
    {
        self.get_array_mut().borrow_mut()
    }
}

impl<C, const N: usize> const AsRef<[C]> for ArrayPolynomial<C, N>
where
    [C; N]: ~const AsRef<[C]>
{
    fn as_ref(&self) -> &[C]
    {
        self.get_array().as_ref()
    }
}
impl<C, const N: usize> const AsMut<[C]> for ArrayPolynomial<C, N>
where
    [C; N]: ~const AsMut<[C]>
{
    fn as_mut(&mut self) -> &mut [C]
    {
        self.get_array_mut().as_mut()
    }
}

impl<C, const N: usize> const ArrayOps<C, N> for ArrayPolynomial<C, N>
{
    type Array<I, const L: usize> = ArrayPolynomial<I, L>;

    fn fill<F>(fill: F) -> Self
    where
        F: ~const FnMut(usize) -> C + ~const Destruct
    {
        <[C; N]>::fill(fill).into()
    }

    fn rfill<F>(fill: F) -> Self
    where
        F: ~const FnMut(usize) -> C + ~const Destruct
    {
        <[C; N]>::rfill(fill).into()
    }

    fn truncate<const M: usize>(self) -> Self::Array<C, M>
    where
        C: ~const Destruct,
        [(); N - M]:
    {
        self.into_array().truncate().into()
    }

    fn rtruncate<const M: usize>(self) -> Self::Array<C, M>
    where
        C: ~const Destruct,
        [(); N - M]:
    {
        self.into_array().rtruncate().into()
    }

    fn resize<const M: usize, F>(self, fill: F) -> Self::Array<C, M>
    where
        F: ~const FnMut(usize) -> C + ~const Destruct,
        C: ~const Destruct
    {
        self.into_array().resize(fill).into()
    }

    fn rresize<const M: usize, F>(self, fill: F) -> Self::Array<C, M>
    where
        F: ~const FnMut(usize) -> C + ~const Destruct,
        C: ~const Destruct
    {
        self.into_array().rresize(fill).into()
    }

    fn extend<const M: usize, F>(self, fill: F) -> Self::Array<C, M>
    where
        F: ~const FnMut(usize) -> C + ~const Destruct,
        [(); M - N]:
    {
        self.into_array().extend(fill).into()
    }

    fn rextend<const M: usize, F>(self, fill: F) -> Self::Array<C, M>
    where
        F: ~const FnMut(usize) -> C + ~const Destruct,
        [(); M - N]:
    {
        self.into_array().rextend(fill).into()
    }

    fn into_const_iter(self) -> array_trait::IntoConstIter<C, N, true>
    {
        self.into_array().into_const_iter()
    }

    fn into_const_iter_reverse(self) -> array_trait::IntoConstIter<C, N, false>
    {
        self.into_array().into_const_iter_reverse()
    }

    fn const_iter(&self) -> array_trait::ConstIter<'_, C, N>
    {
        self.get_array().const_iter()
    }

    fn const_iter_mut(&mut self) -> array_trait::ConstIterMut<'_, C, N>
    {
        self.get_array_mut().const_iter_mut()
    }

    fn map2<M>(self, map: M) -> Self::Array<<M as FnOnce<(C,)>>::Output, N>
    where
        M: ~const FnMut<(C,)> + ~const Destruct
    {
        self.into_array().map2(map).into()
    }

    fn zip2<O, Rhs>(self, other: Rhs) -> Self::Array<(C, O), N>
    where
        Rhs: ~const Into<Self::Array<O, N>>
    {
        self.into_array().zip2(other.into().into_array()).into()
    }

    fn enumerate(self) -> Self::Array<(usize, C), N>
    {
        self.into_array().enumerate().into()
    }

    fn reduce<R>(self, reduce: R) -> Option<C>
    where
        R: ~const FnMut(C, C) -> C + ~const Destruct,
        C: ~const Destruct
    {
        self.into_array().reduce(reduce)
    }

    fn chain<const M: usize, Rhs>(self, rhs: Rhs) -> Self::Array<C, {N + M}>
    where
        Rhs: ~const Into<Self::Array<C, M>>
    {
        self.into_array().chain(rhs.into().into_array()).into()
    }

    fn array_spread<const M: usize>(self) -> ([Self::Array<C, {N / M}>; M], Self::Array<C, {N % M}>)
    where
        [(); M - 1]:,
        C: Copy
    {
        let (spread, right) = self.into_array().array_spread();
        (spread.map2(Into::into), right.into())
    }

    fn array_spread_ref<const M: usize>(&self) -> ([&Self::PaddedArray<C, M, {N / M}>; M], &Self::Array<C, {N % M}>)
    where
        [(); M - 1]:
    {
        let (spread, right) = self.get_array().array_spread_ref();
        (spread.map2(Into::into), right.into())
    }

    fn array_spread_mut<const M: usize>(&mut self) -> ([&mut Self::PaddedArray<C, M, {N / M}>; M], &mut Self::Array<C, {N % M}>)
    where
        [(); M - 1]:
    {
        let (spread, right) = self.get_array_mut().array_spread_mut();
        (spread.map2(Into::into), right.into())
    }

    fn array_rspread<const M: usize>(self) -> (Self::Array<C, {N % M}>, [Self::Array<C, {N / M}>; M])
    where
        [(); M - 1]:,
        C: Copy
    {
        let (left, spread) = self.into_array().array_rspread();
        (left.into(), spread.map2(Into::into))
    }

    fn array_rspread_ref<const M: usize>(&self) -> (&Self::Array<C, {N % M}>, [&Self::PaddedArray<C, M, {N / M}>; M])
    where
        [(); M - 1]:
    {
        let (left, spread) = self.get_array().array_rspread_ref();
        (left.into(), spread.map2(Into::into))
    }

    fn array_rspread_mut<const M: usize>(&mut self) -> (&mut Self::Array<C, {N % M}>, [&mut Self::PaddedArray<C, M, {N / M}>; M])
    where
        [(); M - 1]:
    {
        let (left, spread) = self.get_array_mut().array_rspread_mut();
        (left.into(), spread.map2(Into::into))
    }

    fn array_spread_exact<const M: usize>(self) -> [Self::Array<C, {N / M}>; M]
    where
        [(); M - 1]:,
        [(); 0 - N % M]:
    {
        self.into_array().array_spread_exact().map2(Into::into)
    }

    fn array_spread_exact_ref<const M: usize>(&self) -> [&Self::PaddedArray<C, M, {N / M}>; M]
    where
        [(); M - 1]:,
        [(); 0 - N % M]:
    {
        self.get_array().array_spread_exact_ref().map2(Into::into)
    }

    fn array_spread_exact_mut<const M: usize>(&mut self) -> [&mut Self::PaddedArray<C, M, {N / M}>; M]
    where
        [(); M - 1]:,
        [(); 0 - N % M]:
    {
        self.get_array_mut().array_spread_exact_mut().map2(Into::into)
    }

    fn array_chunks<const M: usize>(self) -> ([Self::Array<C, M>; N / M], Self::Array<C, {N % M}>)
    {
        RCurry::<fn([C; M]) -> Self::Array<C, M>, ([[C; M]; N / M],)>::rcurry(
            ArrayOps::map2,
            From::from
        ).fn_zip(From::from).call_once(self.into_array().array_chunks())
    }

    fn array_chunks_ref<const M: usize>(&self) -> ([&Self::Array<C, M>; N / M], &Self::Array<C, {N % M}>)
    {
        let (chunks, left) = self.get_array().array_chunks_ref();
        (chunks.map2(Into::into), left.into())
    }

    fn array_chunks_mut<const M: usize>(&mut self) -> ([&mut Self::Array<C, M>; N / M], &mut Self::Array<C, {N % M}>)
    {
        let (chunks, left) = self.get_array_mut().array_chunks_mut();
        (chunks.map2(Into::into), left.into())
    }

    fn array_rchunks<const M: usize>(self) -> (Self::Array<C, {N % M}>, [Self::Array<C, M>; N / M])
    {
        Into::into.fn_zip(
            RCurry::<fn([C; M]) -> Self::Array<C, M>, ([[C; M]; N / M],)>::rcurry(
                ArrayOps::map2,
                Into::into
            )
        ).call(self.into_array().array_rchunks())
    }

    fn array_rchunks_ref<const M: usize>(&self) -> (&Self::Array<C, {N % M}>, [&Self::Array<C, M>; N / M])
    {
        let (right, chunks) = self.get_array().array_rchunks_ref();
        (right.into(), chunks.map2(Into::into))
    }

    fn array_rchunks_mut<const M: usize>(&mut self) -> (&mut Self::Array<C, {N % M}>, [&mut Self::Array<C, M>; N / M])
    {
        let (right, chunks) = self.get_array_mut().array_rchunks_mut();
        (right.into(), chunks.map2(Into::into))
    }

    fn array_chunks_exact<const M: usize>(self) -> [Self::Array<C, M>; N / M]
    where
        [(); 0 - N % M]:,
        [(); N / M]:
    {
        self.into_array().array_chunks_exact().map2(Into::into)
    }

    fn array_chunks_exact_ref<const M: usize>(&self) -> [&Self::Array<C, M>; N / M]
    where
        [(); 0 - N % M]:,
        [(); N / M]:
    {
        self.get_array().array_chunks_exact_ref().map2(Into::into)
    }

    fn array_chunks_exact_mut<const M: usize>(&mut self) -> [&mut Self::Array<C, M>; N / M]
    where
        [(); 0 - N % M]:,
        [(); N / M]:
    {
        self.get_array_mut().array_chunks_exact_mut().map2(Into::into)
    }

    fn split_array<const M: usize>(self) -> (Self::Array<C, M>, Self::Array<C, {N - M}>)
    where
        [(); N - M]:
    {
        Into::into.fn_zip(Into::into)
            .call(self.into_array().split_array())
    }

    fn split_array_ref2<const M: usize>(&self) -> (&Self::Array<C, M>, &Self::Array<C, {N - M}>)
    where
        [(); N - M]:
    {
        Into::into.fn_zip(Into::into)
            .call(self.get_array().split_array_ref2())
    }

    fn split_array_mut2<const M: usize>(&mut self) -> (&mut Self::Array<C, M>, &mut Self::Array<C, {N - M}>)
    where
        [(); N - M]:
    {
        Into::into.fn_zip(Into::into)
            .call(self.get_array_mut().split_array_mut2())
    }

    fn rsplit_array<const M: usize>(self) -> (Self::Array<C, {N - M}>, Self::Array<C, M>)
    where
        [(); N - M]:
    {
        Into::into.fn_zip(Into::into)
            .call(self.into_array().rsplit_array())
    }

    fn rsplit_array_ref2<const M: usize>(&self) -> (&Self::Array<C, {N - M}>, &Self::Array<C, M>)
    where
        [(); N - M]:
    {
        Into::into.fn_zip(Into::into)
            .call(self.get_array().rsplit_array_ref2())
    }

    fn rsplit_array_mut2<const M: usize>(&mut self) -> (&mut Self::Array<C, {N - M}>, &mut Self::Array<C, M>)
    where
        [(); N - M]:
    {
        Into::into.fn_zip(Into::into)
            .call(self.get_array_mut().rsplit_array_mut2())
    }
}