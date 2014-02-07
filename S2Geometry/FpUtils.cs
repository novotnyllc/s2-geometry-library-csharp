using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Google.Common.Geometry
{
    // Data from here
    //http://grepcode.com/file_/repository.grepcode.com/java/root/jdk/openjdk/6-b14/sun/misc/DoubleConsts.java/?v=source
    class DoubleConsts
    {
         /**
     * The number of logical bits in the significand of a
     * <code>double</code> number, including the implicit bit.
     */
    public const int SIGNIFICAND_WIDTH   = 53;

    public static readonly double MIN_NORMAL = BitConverter.Int64BitsToDouble(0x0010000000000000L);

       /**
     * Maximum exponent a finite <code>double</code> number may have.
     * It is equal to the value returned by
     * <code>Math.ilogb(Double.MAX_VALUE)</code>.
     */
    public const int     MAX_EXPONENT    = 1023;

    /**
     * Minimum exponent a normalized <code>double</code> number may
     * have.  It is equal to the value returned by
     * <code>Math.ilogb(Double.MIN_NORMAL)</code>.
     */
    public const int     MIN_EXPONENT    = -1022;

            /**
     * Bit mask to isolate the exponent field of a
     * <code>double</code>.
     */
    public const long    EXP_BIT_MASK    = 0x7FF0000000000000L;

         /**
     * Bias used in representing a <code>double</code> exponent.
     */
    public const int     EXP_BIAS        = 1023;
    }

    // Data from here http://grepcode.com/file_/repository.grepcode.com/java/root/jdk/openjdk/6-b14/sun/misc/FpUtils.java/?v=source
    class FpUtils
    {
        static readonly double twoToTheDoubleScaleUp = powerOfTwoD(512);
        static readonly double twoToTheDoubleScaleDown = powerOfTwoD(-512);

        /**
    * Returns a floating-point power of two in the normal range.
    */
        static double powerOfTwoD(int n)
        {
            Debug.Assert(n >= DoubleConsts.MIN_EXPONENT && n <= DoubleConsts.MAX_EXPONENT);
            return BitConverter.Int64BitsToDouble((((long)n + (long)DoubleConsts.EXP_BIAS) <<
                                            (DoubleConsts.SIGNIFICAND_WIDTH - 1))
                                           & DoubleConsts.EXP_BIT_MASK);
        }
        /**
    * Return <code>d</code> &times;
    * 2<sup><code>scale_factor</code></sup> rounded as if performed
    * by a single correctly rounded floating-point multiply to a
    * member of the double value set.  See <a
    * href="http://java.sun.com/docs/books/jls/second_edition/html/typesValues.doc.html#9208">&sect;4.2.3</a>
    * of the <a href="http://java.sun.com/docs/books/jls/html/">Java
    * Language Specification</a> for a discussion of floating-point
    * value sets.  If the exponent of the result is between the
    * <code>double</code>'s minimum exponent and maximum exponent,
    * the answer is calculated exactly.  If the exponent of the
    * result would be larger than <code>doubles</code>'s maximum
    * exponent, an infinity is returned.  Note that if the result is
    * subnormal, precision may be lost; that is, when <code>scalb(x,
    * n)</code> is subnormal, <code>scalb(scalb(x, n), -n)</code> may
    * not equal <i>x</i>.  When the result is non-NaN, the result has
    * the same sign as <code>d</code>.
    *
    *<p>
    * Special cases:
    * <ul>
    * <li> If the first argument is NaN, NaN is returned.
    * <li> If the first argument is infinite, then an infinity of the
    * same sign is returned.
    * <li> If the first argument is zero, then a zero of the same
    * sign is returned.
    * </ul>
    *
    * @param d number to be scaled by a power of two.
    * @param scale_factor power of 2 used to scale <code>d</code>
    * @return <code>d * </code>2<sup><code>scale_factor</code></sup>
    * @author Joseph D. Darcy
    */
        public static double scalb(double d, int scale_factor) {
        /*
         * This method does not need to be declared strictfp to
         * compute the same correct result on all platforms.  When
         * scaling up, it does not matter what order the
         * multiply-store operations are done; the result will be
         * finite or overflow regardless of the operation ordering.
         * However, to get the correct result when scaling down, a
         * particular ordering must be used.
         *
         * When scaling down, the multiply-store operations are
         * sequenced so that it is not possible for two consecutive
         * multiply-stores to return subnormal results.  If one
         * multiply-store result is subnormal, the next multiply will
         * round it away to zero.  This is done by first multiplying
         * by 2 ^ (scale_factor % n) and then multiplying several
         * times by by 2^n as needed where n is the exponent of number
         * that is a covenient power of two.  In this way, at most one
         * real rounding error occurs.  If the double value set is
         * being used exclusively, the rounding will occur on a
         * multiply.  If the double-extended-exponent value set is
         * being used, the products will (perhaps) be exact but the
         * stores to d are guaranteed to round to the double value
         * set.
         *
         * It is _not_ a valid implementation to first multiply d by
         * 2^MIN_EXPONENT and then by 2 ^ (scale_factor %
         * MIN_EXPONENT) since even in a strictfp program double
         * rounding on underflow could occur; e.g. if the scale_factor
         * argument was (MIN_EXPONENT - n) and the exponent of d was a
         * little less than -(MIN_EXPONENT - n), meaning the final
         * result would be subnormal.
         *
         * Since exact reproducibility of this method can be achieved
         * without any undue performance burden, there is no
         * compelling reason to allow double rounding on underflow in
         * scalb.
         */

        // magnitude of a power of two so large that scaling a finite
        // nonzero value by it would be guaranteed to over or
        // underflow; due to rounding, scaling down takes takes an
        // additional power of two which is reflected here
        int MAX_SCALE = DoubleConsts.MAX_EXPONENT + -DoubleConsts.MIN_EXPONENT +
                              DoubleConsts.SIGNIFICAND_WIDTH + 1;
        int exp_adjust = 0;
        int scale_increment = 0;
        double exp_delta = Double.NaN;

        // Make sure scaling factor is in a reasonable range

        if(scale_factor < 0) {
            scale_factor = Math.Max(scale_factor, -MAX_SCALE);
            scale_increment = -512;
            exp_delta = twoToTheDoubleScaleDown;
        }
        else {
            scale_factor = Math.Min(scale_factor, MAX_SCALE);
            scale_increment = 512;
            exp_delta = twoToTheDoubleScaleUp;
        }

        // Calculate (scale_factor % +/-512), 512 = 2^9, using
        // technique from "Hacker's Delight" section 10-2.
        uint u = unchecked ((uint)(scale_factor >> 9 - 1));
        u =   u >> 32 - 9;
            int t = unchecked ((int)u);
        exp_adjust = ((scale_factor + t) & (512 -1)) - t;

        d *= powerOfTwoD(exp_adjust);
        scale_factor -= exp_adjust;

        while(scale_factor != 0) {
            d *= exp_delta;
            scale_factor -= scale_increment;
        }
        return d;
    }


    }
}
