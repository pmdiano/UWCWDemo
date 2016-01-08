using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace DemoUWCW
{
    class InterfaceUWCW
    {
        [DllImport("CoreUWCW.dll", CallingConvention = CallingConvention.Cdecl)]
        private unsafe static extern short detectUWCW(float* RP_dBm,
        short M, short N, float* Ang_plot, float Ang_endpt, float fft_dR,
        float MagMax, float MagMin, bool isFirstFrame,
        short* rangeIndices, short* angleIndices);

        [DllImport("CoreUWCW.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void InitUWCW();

        [DllImport("CoreUWCW.dll", CallingConvention = CallingConvention.Cdecl)]
        private unsafe static extern void getScore(float* currentScore, float* prevScore,
        float* minRange, float* maxRange);

        [DllImport("CoreUWCW.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void SetAlgorithm(bool useOriginal);


        public static short MainInterface(float[] RP_dBm, short M, short N, float[] Ang_plot,
            float Ang_endpt, float fft_dR, float MagMax, float MagMin, bool isFirstFrame,
            short[] ranges, short[] angles)
        {
            short totalPixNum;

            unsafe
            {
                fixed (float* fptr = RP_dBm, aptr = Ang_plot)
                {
                    fixed (short* rangePtr = ranges, anglePtr = angles)
                    {
                        totalPixNum = detectUWCW(fptr, M, N, aptr, Ang_endpt, fft_dR, MagMax, MagMin,
                            isFirstFrame, rangePtr, anglePtr);
                    }
                }
            }

            return totalPixNum;
        }

        public static void GetScoreRange(ref float score, ref float prevScore, ref float minRange, ref float maxRange)
        {
            unsafe
            {
                fixed (float* scorePtr = &score, prevScorePtr = &prevScore)
                {
                    fixed (float* minRPtr = &minRange, maxRPtr = &maxRange)
                    {
                        getScore(scorePtr, prevScorePtr, minRPtr, maxRPtr);
                    }
                }
            }
        }
    }
}