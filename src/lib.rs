#![cfg_attr(not(any(test, feature = "std")), no_std)]

#![feature(const_trait_impl)]
#![feature(const_mut_refs)]
#![feature(allocator_api)]

#![feature(generic_const_exprs)]
#![feature(const_closures)]

moddef::moddef!(
    flat(pub) mod {
        chebyshev_polynomial,
        mul_polynomial,
        plot for cfg(test),
        polynomial_nd,
        polynomial,
        product_polynomial
    }
);

#[cfg(test)]
mod tests {
    use array__ops::ArrayOps;

    use super::*;

    #[cfg(feature = "std")]
    #[test]
    fn mul()
    {
        let p = vec![1.0, -1.0];
        //println!("{:?}", p.clone().mul_polynomial(p.clone()));

        const P: [f32; 2] = [1.0, 1.0];
        let q: [f32; 3] = P.mul_polynomial(P);
        println!("{:?}", q);

        let q = p.as_slice().mul_polynomial(P);
        println!("{:?}", q);
    }

    #[cfg(feature = "std")]
    #[test]
    fn eval_nd()
    {
        use array__ops::ArrayNdOps;

        const A: [[u128; 3]; 3] = [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 1]
        ];

        let xy: [[[u128; 2]; 3]; 3] = ArrayNdOps::fill_nd(|i| i.map2(const |i| i as u128));

        assert_eq!(xy, [
            [[0, 0], [0, 1], [0, 2]],
            [[1, 0], [1, 1], [1, 2]],
            [[2, 0], [2, 1], [2, 2]]
        ]);

        println!("{:?}", xy.map_nd(|xy: [u128; 2]| A.evaluate_as_polynomial_nd(xy)))
    }

    #[cfg(feature = "std")]
    #[test]
    fn it_works()
    {
        type T = i128;
        const KIND: usize = 2;
        const ORDER: usize = 13;
        const N: usize = ORDER + 1;

        println!("KIND = {}\nORDER = {}", KIND, ORDER);

        let cheb_vec: Vec<T> = ChebyshevPolynomial::new(KIND, ORDER).into();
        println!("cheb = {:?}", cheb_vec);
        let cheb: [T; N] = ChebyshevPolynomial::new(KIND, ORDER).try_into().ok().unwrap();
        println!("CHEB = {:?}", cheb);

        let max: T = cheb.reduce(Ord::max).unwrap();
        let min: T = cheb.reduce(Ord::min).unwrap();

        print!("max = {}", max);
        if max > i64::MAX as i128 {print!(" > i64")} // N < 54
        if max > i32::MAX as i128 {print!(" > i32")} // N < 28
        if max > i16::MAX as i128 {print!(" > i16")} // N < 14
        if max > i8::MAX as i128 {print!(" > i8")} // N < 9
        println!("");
        print!("min = {}", min);
        if min < i64::MIN as i128 {print!(" < i64")}
        if min < i32::MIN as i128 {print!(" < i32")}
        if min < i16::MIN as i128 {print!(" < i16")}
        if min < i8::MIN as i128 {print!(" < i8")}
        println!("");
        
        assert_eq!(cheb.to_vec(), cheb_vec);

        const TEST_COUNT: usize = (i8::MAX as i64 - i8::MIN as i64 + 1) as usize;
        let x: [T; TEST_COUNT] = ArrayOps::fill(const |i| i as T + i8::MIN as T);

        let y: [T; TEST_COUNT] = x.map(|x| cheb.evaluate_as_polynomial(x));
        let y_vec = x.map(|x| cheb_vec.clone().evaluate_as_polynomial(x));

        assert_eq!(y, y_vec);
    }

    mod plot
    {
        use crate::plot::*;

        type T = f32;
        
        const PLOT_TARGET: &str = "plots";

        use crate::{ChebyshevPolynomial, Polynomial, PolynomialNd};
        use array__ops::ArrayOps;
        use currying::Curry;
        use linspace::LinspaceArray;

        #[test]
        fn plot_cheb()
        {
            const KIND: usize = 2;
            const ORDER: usize = 2;
            const N: usize = ORDER + 1;
            
            let cheb: [T; N] = ChebyshevPolynomial::new(KIND, ORDER).try_into().ok().unwrap();
            
            const X0: T = 0.0;
            const X1: T = 4.0;

            const RESOLUTION: usize = 1024;
            let x: [T; RESOLUTION] = (X0..X1).linspace_array();

            let y: [T; RESOLUTION] = x.map2(Polynomial::evaluate_as_polynomial.curry(cheb));
            
            let plot_title: &str = &format!("{ORDER}. order Chebyshev of the {KIND}. kind, x = {X0}..{X1}:");
            let plot_path: &str = &format!("{PLOT_TARGET}/chebyshev.png"); //&format!("{PLOT_TARGET}/chebyshev_k{KIND}_o{ORDER}_x{X0}_{X1}.png");

            plot_curve(plot_title, plot_path, x, y).unwrap()
        }
        
        #[test]
        fn plot_polynomial()
        {
            //const POLYNOMIAL: [T; 5] = [1.0, 0.0, -10.0, 0.0, 1.0];
            //const POLYNOMIAL: [T; 6] = [2.0, -4.0, 0.0, 0.0, 0.0, 1.0];
            //const POLYNOMIAL: [T; 3] = [1.0, 5.0, 6.0];
            const POLYNOMIAL: [T; 6] = [1.0, -1.0, 0.0, 0.0, 0.0, 1.0];
            const ORDER: usize = POLYNOMIAL.len() + 1;
            
            const X0: T = -5 as T;
            const X1: T = 5 as T;

            const RESOLUTION: usize = 1024;
            let x: [T; RESOLUTION] = (X0..X1).linspace_array();

            let y: [T; RESOLUTION] = x.map(Polynomial::evaluate_as_polynomial.curry(POLYNOMIAL));
            
            let plot_title: &str = &format!("{ORDER}. order polynomial, x = {X0}..{X1}:");
            let plot_path: &str = &format!("{PLOT_TARGET}/polynomial.png"); //&format!("{PLOT_TARGET}/polynomial_{:?}_x{X0}_{X1}.png", POLYNOMIAL);

            plot_curve(plot_title, plot_path, x, y).unwrap()
        }
        #[test]
        fn plot_polynomial_2d_parametric()
        {
            const POLYNOMIAL: [[[T; 2]; 2]; 3] = [
                [
                    [0., 1.],
                    [0., 1.]
                ],
                [
                    [0., -1.],
                    [1., -3.]
                ],
                [
                    [0., 1.],
                    [1., 0.]
                ]
            ];
            
            const X0: T = -10 as T;
            const X1: T = 10 as T;

            const Y0: T = -10 as T;
            const Y1: T = 10 as T;

            const RESOLUTION: [usize; 2] = [32, 32];
            let x: [T; RESOLUTION[0]] = (X0..X1).linspace_array();
            let y: [T; RESOLUTION[1]] = (Y0..Y1).linspace_array();
            
            let plot_title: &str = &format!("2D polynomial, (x, y) = ({X0}..{X1}, {Y0}..{Y1}):");
            let plot_path: &str = &format!("{PLOT_TARGET}/polynomial_2d_parametric.svg");

            plot_parametric_curve_2d(
                plot_title, plot_path, x, y,
                |x, y| POLYNOMIAL.map2(|p| p.evaluate_as_polynomial_nd([x, y]))
            ).unwrap()
        }
        
        #[test]
        fn plot_polynomial_2d()
        {
            const POLYNOMIAL: [[T; 3]; 3] = [
                [2., 1., -0.],
                [0., 0., 1.],
                [-1., 0., -1.]
            ];
            
            const X0: T = -10 as T;
            const X1: T = 10 as T;

            const Y0: T = -10 as T;
            const Y1: T = 10 as T;

            const RESOLUTION: [usize; 2] = [32, 32];
            let x: [T; RESOLUTION[0]] = (X0..X1).linspace_array();
            let y: [T; RESOLUTION[1]] = (Y0..Y1).linspace_array();
            
            let plot_title: &str = &format!("2D polynomial, (x, y) = ({X0}..{X1}, {Y0}..{Y1}):");
            let plot_path: &str = &format!("{PLOT_TARGET}/polynomial_2d.svg"); //&format!("{PLOT_TARGET}/polynomial_{:?}_x{X0}_{X1}.png", POLYNOMIAL);

            plot_curve_2d(
                plot_title, plot_path, x, y,
                |x, y| POLYNOMIAL.evaluate_as_polynomial_nd([x, y])
            ).unwrap()
        }

        #[test]
        fn plot_polynomial_rad()
        {
            //const POLYNOMIAL: [T; 5] = [1.0, 0.0, -10.0, 0.0, 1.0];
            const POLYNOMIAL: [[T; 5]; 2] = [
                [1., 0., -10., 0., 1.],
                [1000000., 10., 0., 0., 0.]
            ];

            const R0: T = 0.0;
            const R1: T = 1.0;

            const THETA0: T = 0.0;
            const THETA1: T = core::f32::consts::TAU*4.0;

            const RESOLUTION: [usize; 2] = [16, 128];
            let r: [T; RESOLUTION[0]] = (R0..R1).linspace_array();
            let theta: [T; RESOLUTION[1]] = (THETA0..THETA1).linspace_array();
            
            let plot_title: &str = &format!("2D radial polynomial, (r, theta) = ({R0}..{R1}, {THETA0}..{THETA1}):");
            let plot_path: &str = &format!("{PLOT_TARGET}/polynomial_2d_rad.svg");

            plot_curve_2d_rad(plot_title, plot_path, r, theta, |r, theta| POLYNOMIAL.evaluate_as_polynomial_nd([r, theta])).unwrap()
        }
    
        #[test]
        fn plot_polynomial_parametric_rad()
        {
            const POLYNOMIAL: [[[T; 2]; 2]; 3] = [
                [
                    [0., 1.],
                    [0.1, 0.]
                ],
                [
                    [0., -1.5],
                    [1., 0.1]
                ],
                [
                    [0., 0.4],
                    [0.2, 1.]
                ]
            ];
    
            const R0: T = 2.0;
            const R1: T = 5.0;
    
            const THETA0: T = 0.0;
            const THETA1: T = core::f32::consts::TAU;
    
            const RESOLUTION: [usize; 2] = [16, 128];
            let r: [T; RESOLUTION[0]] = (R0..R1).linspace_array();
            let theta: [T; RESOLUTION[1]] = (THETA0..THETA1).linspace_array();
            
            let plot_title: &str = &format!("2D radial polynomial, (r, theta) = ({R0}..{R1}, {THETA0}..{THETA1}):");
            let plot_path: &str = &format!("{PLOT_TARGET}/polynomial_2d_parametric_rad.svg");
    
            plot_parametric_curve_2d_rad(
                plot_title, plot_path, r, theta,
                |r, theta| POLYNOMIAL.map2(|p| p.evaluate_as_polynomial_nd([r, theta]))
            ).unwrap()
        }
    }

    #[cfg(disabled)]
    mod size
    {
        use num::Float;

        #[test]
        fn measure_size()
        {
            const N: usize = 64;
            println!("[");
            for order in 0..N
            {
                println!("[{}, {}, {}]", order, max_bits_of_cheb_coeff(1, order), guess(1, order));
            }
            println!("]")
        }
    
        fn guess(kind: usize, order: usize) -> u32
        {
            (
                2.0.powi(2*order as i32 - kind as i32 - 2)
            ).abs().log2().ceil() as u32
            //ceil(log2(|(-1)^(M+N) * 2^(2N-M-1)|))
            //ceil(log2(|2^(2N-M-1)|))
        }
    
        fn max_bits_of_cheb_coeff(kind: usize, order: usize) -> u32
        {
            let cheb = Vec::<i128>::new_chebyshev_polynomial_of_kind(kind, order).unwrap();
            
            let max = cheb.clone().into_iter().reduce(Ord::max).unwrap();
            let min = cheb.into_iter().reduce(Ord::min).unwrap();
            
            let mut bits = 1;
            let mut n = max.abs().max(min.abs());
            while n != 0
            {
                n >>= 1;
                bits += 1;
            }
            bits
        }
    }
}