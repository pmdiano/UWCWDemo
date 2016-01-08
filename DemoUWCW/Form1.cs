using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;
using System.Media;
using System.Diagnostics;

using csmatio.types;
using csmatio.io;

namespace DemoUWCW
{
    public partial class Form1 : Form
    {
        bool justShowOriginal = false;
        bool useOriginalAlgrm = false;

        // reading file variables
        private string firstName;
        private string parameterName;
        private string currentName;
        private string frameType;
        private int firstId;
        private int currentId;
        private int lastId;
        private bool loaded;
        private StreamWriter file;

        // Prm_ parameters
        float MagMax;
        float MagMin;
        int fft_bins;        
        float fft_dR;
        float Ang_endpt;
        Int16 Az_bins = 0;
        Int16[] Az_cnt = null;
        Single[] RP_dBm = null;
        Single[] Ang_plot = null;

        // output variables
        float score = 0;
        float prevScore = 0;
        float minRange = 0;
        float maxRange = 0;

        // image display parameters
        byte[,] COLORMAP_JET;
        short[] ranges = null;
        short[] angles = null;
        short pixNum = 0;
        bool cablesShown = false;

        // timing functions
        int totalFrames = 0;
        double totalTime = 0.0;
        Stopwatch sw = null;

        public Form1()
        {
            InitializeComponent();
            this.FormBorderStyle = FormBorderStyle.Fixed3D;
            initColormap();
            loaded = false;
        }

        // initialize the jet colormap
        private void initColormap()
        {
            COLORMAP_JET = new byte[64, 3]{
                {0,     0,   143},
                {0,     0,   159},
                {0,     0,   175},
                {0,     0,   191},
                {0,     0,   207},
                {0,     0,   223},
                {0,     0,   239},
                {0,     0,   255},
                {0,    16,   255},
                {0,    32,   255},
                {0,    48,   255},
                {0,    64,   255},
                {0,    80,   255},
                {0,    96,   255},
                {0,   112,   255},
                {0,   128,   255},
                {0,   143,   255},
                {0,   159,   255},
                {0,   175,   255},
                {0,   191,   255},
                {0,   207,   255},
                {0,   223,   255},
                {0,   239,   255},
                {0,   255,   255},
                {16,   255,   239},
                {32,   255,   223},
                {48,   255,   207},
                {64,   255,   191},
                {80,   255,   175},
                {96,   255,   159},
                {112,   255,   143},
                {128,   255,   128},
                {143,   255,   112},
                {159,   255,    96},
                {175,   255,    80},
                {191,   255,    64},
                {207,   255,    48},
                {223,   255,    32},
                {239,   255,    16},
                {255,   255,     0},
                {255,   239,     0},
                {255,   223,     0},
                {255,   207,     0},
                {255,   191,     0},
                {255,   175,     0},
                {255,   159,     0},
                {255,   143,     0},
                {255,   128,     0},
                {255,   112,     0},
                {255,    96,     0},
                {255,    80,     0},
                {255,    64,     0},
                {255,    48,     0},
                {255,    32,     0},
                {255,    16,     0},
                {255,     0,     0},
                {239,     0,     0},
                {223,     0,     0},
                {207,     0,     0},
                {191,     0,     0},
                {175,     0,     0},
                {159,     0,     0},
                {143,     0,     0},
                {128,     0,     0}
            };
        }

        private void btnStep_Click(object sender, EventArgs e)
        {
            if (!loaded)
            {
                MessageBox.Show("Please load data first.");
                return;
            }
            try
            {
                processOneFrame(sender, e);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Next Frame error: " + ex.Message);
            }
        }

        private void btnRun_Click(object sender, EventArgs e)
        {
            if (!loaded)
            {
                MessageBox.Show("Please load data first.");
                return;
            }
            try
            {
                Application.Idle += processOneFrame;
            }
            catch (Exception ex)
            {
                MessageBox.Show("All Frames error: " + ex.Message);
            }
        }

        private void btnPause_Click(object sender, EventArgs e)
        {
            try
            {
                Application.Idle -= processOneFrame;
            }
            catch (Exception ex)
            {
                MessageBox.Show("Pause error: " + ex.Message);
            }
        }

        private void btnExit_Click(object sender, EventArgs e)
        {
            this.Close();
        }

        private void exitToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.Close();
        }

        private void openToolStripMenuItem_Click(object sender, EventArgs e)
        {
            OpenFileDialog firstFrameDialog = new OpenFileDialog();
            firstFrameDialog.Title = "Open a frame";
            firstFrameDialog.Filter = "mat frame (*.mat)|*.mat|bfr frame (*.bfr)|*.bfr";
            firstFrameDialog.FilterIndex = 2;
            firstFrameDialog.RestoreDirectory = true;
            if (!loaded)
            {
                firstFrameDialog.InitialDirectory = Directory.GetCurrentDirectory();
            }

            if (firstFrameDialog.ShowDialog() == DialogResult.OK)
            {
                ranges = null;
                angles = null;
                pixNum = 0;

                try
                {
                    // get first name, file type, and paramter name
                    firstName = firstFrameDialog.FileName;                    
                    frameType = firstName.Substring(firstName.Length - 3);

                    // get first and current ID
                    string str = firstName.Substring(firstName.Length - 8, 4);
                    firstId = Int32.Parse(str);
                    currentId = firstId;

                    // get last ID
                    string[] fileNames = Directory.GetFiles(Path.GetDirectoryName(firstName), "*." + frameType);
                    int i = fileNames.Length - 1;
                    while (fileNames[i].Length != firstName.Length || fileNames[i].Substring(0, firstName.Length-8) != firstName.Substring(0, firstName.Length-8))
                    {
                        i--;
                    }
                    str = fileNames[i].Substring(firstName.Length - 8, 4);
                    lastId = Int32.Parse(str);

                    // get result file stream writer
                    if (file != null) file.Close();
                    file = new StreamWriter("frameScore.txt");

                    // reset timer
                    totalFrames = 0;
                    totalTime = 0.0;
                    sw = new Stopwatch();

                    // show first frame
                    readParameters();
                    InterfaceUWCW.InitUWCW();
                    resetVariables();                    
                    processOneFrame(sender, e);
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Error in openning the first frame: " + ex.Message);
                }

                loaded = true;
                diaplayRangeParameters();
                groupBox2.Text = "Dataset\n" + firstName.Substring(0, firstName.Length - 8);                
            }           
        }

        private void resetVariables()
        {
            score = 0;
            prevScore = 0;
            minRange = 0;
            maxRange = 0;
        }

        private void diaplayRangeParameters()
        {
            labelStart.Text = "0 m";
            labelEnd.Text = string.Format("{0:0.00} m", fft_bins * fft_dR);

            labelAngleLeft.BorderStyle = BorderStyle.Fixed3D;
            labelAngleMiddle.BorderStyle = BorderStyle.Fixed3D;
            labelAngleRight.BorderStyle = BorderStyle.Fixed3D;

            toolStripStatusLabel2.Text = "MagMin: " + MagMin.ToString();
            toolStripStatusLabel3.Text = "MagMax: " + MagMax.ToString();
            toolStripStatusLabel4.Text = "fft_dR: " + string.Format("{0:0.000}", fft_dR) + " m";
        }

        // display angle parameters, starting and ending positions;
        // also display size information
        private void displayAngleParameters(short[] Az_cnt, float[] Ang_plot)
        {
            if (Ang_plot != null)
            {
                toolStripStatusLabel1.Text = "Size: " + fft_bins.ToString() + " X " + Ang_plot.Length.ToString();
                labelAngleMiddle.Text = Ang_plot[Ang_plot.Length / 2].ToString();
                string first = string.Format("{0:0.00}", Ang_plot[0]);
                string last = string.Format("{0:0.00}", Ang_plot[Ang_plot.Length - 1]);
                if (currentId % 2 != 0)
                {
                    labelAngleLeft.Text = first;
                    labelAngleRight.Text = last;
                }
                else
                {
                    labelAngleLeft.Text = last;
                    labelAngleRight.Text = first;
                }
            }
            else if (Az_cnt != null)
            {
                toolStripStatusLabel1.Text = "Size: " + fft_bins.ToString() + " X " + Az_cnt.Length.ToString();
                labelAngleMiddle.Text = Az_cnt[Az_cnt.Length / 2].ToString();
                string first = Az_cnt[0].ToString();
                string last = Az_cnt[Az_cnt.Length - 1].ToString();
                if (currentId % 2 != 0)
                {
                    labelAngleLeft.Text = first;
                    labelAngleRight.Text = last;
                }
                else
                {
                    labelAngleLeft.Text = last;
                    labelAngleRight.Text = first;
                }
            }            
            
            toolStripStatusLabel5.Text = "Frame #: " + currentId.ToString();
        }

        // set native data parameters: MagMin, MagMax, etc.
        private void readParameters()
        {
            string pathName = null, pureFirstName = null;
            for (int i = firstName.Length - 1; i >= 0; i--)
            {
                if (firstName[i] == '\\')
                {
                    pathName = firstName.Substring(0, i + 1);
                    pureFirstName = firstName.Substring(i + 1);
                    break;
                }
            }

            parameterName = pathName + "Prm_" + pureFirstName;
            parameterName = parameterName.Substring(0, parameterName.Length - 8);

            switch (frameType)
            {
                case "bfr":
                    // for .bfr data, the name as an additional "_" which we need to remove
                    parameterName = parameterName.Substring(0, parameterName.Length - 1);
                    parameterName += ".m";
                    readParametersBfr(parameterName, out MagMax, out MagMin, out fft_bins, out fft_dR);
                    break;
                case "mat":
                    parameterName += ".mat";
                    readParametersMat(parameterName, out Ang_endpt, out MagMax, out MagMin, out fft_dR);
                    break;
                default:
                    break;
            }
        }

        private void readParametersBfr(string parameterName, out float MagMax, out float MagMin, out int fft_bins, out float fft_dR)
        {
            MagMax = 0;
            MagMin = 0;
            fft_dR = 0;
            fft_bins = 0;
            string[] stringSep = new string[] {"=", ";", " ", "\t"};

            string line;
            StreamReader file = new StreamReader(parameterName);
            while ((line = file.ReadLine()) != null)
            {
                string[] words = line.Split(stringSep, StringSplitOptions.RemoveEmptyEntries);
                if (words.Length == 0) continue;
                switch (words[0])
                {
                    case "MagMax":
                        MagMax = float.Parse(words[1]);
                        break;
                    case "MagMin":
                        MagMin = float.Parse(words[1]);
                        break;
                    case "fft_bins":
                        fft_bins = int.Parse(words[1]);
                        break;
                    case "fft_dR":
                        fft_dR = float.Parse(words[1]);
                        break;
                    default:
                        break;
                }
            }
            file.Close();
        }

        private void readParametersMat(string parameterName, out float Ang_endpt, out float MagMax, out float MagMin, out float fft_dR)
        {
            Ang_endpt = 0;
            MagMax = 0;
            MagMin = 0;
            fft_dR = 0;

            MatFileReader mfr = new MatFileReader(parameterName);
            foreach (MLArray mla in mfr.Data)
            {
                if (mla is MLDouble)
                {
                    MLDouble mld = mla as MLDouble;
                    string name = System.Text.ASCIIEncoding.ASCII.GetString(mld.GetNameToByteArray());
                    switch (name)
                    {
                        case "Ang_endpt":
                            Ang_endpt = (float)mld.Get(0, 0);
                            break;
                        case "MagMax":
                            MagMax = (float)mld.Get(0, 0);
                            break;
                        case "MagMin":
                            MagMin = (float)mld.Get(0, 0);
                            break;
                        case "fft_dR":
                            fft_dR = (float)mld.Get(0, 0);
                            break;
                        default:
                            break;
                    }
                }
            }            
        }

        // open one original image and show it
        private void processOneFrame(object sender, EventArgs e)
        {
            if (currentId > lastId)
            {
                MessageBox.Show("reaches the last frame");
                Application.Idle -= processOneFrame;                
                return;
            }

            // get name and read in one frame
            currentName = firstName.Substring(0, firstName.Length - 8) + currentId.ToString().PadLeft(4, '0') + "." + frameType;
            switch (frameType)
            {
                case "bfr":
                    readOneFrameBfr(currentName, out Az_bins, out Az_cnt, out RP_dBm);
                    Ang_plot = cnt2Ang(Az_cnt);
                    Ang_endpt = Math.Max(Math.Abs(Ang_plot[0]), Math.Abs(Ang_plot[Ang_plot.Length - 1]));
                    RP_dBm = flipColToRow(RP_dBm, (int)Az_bins, (int)fft_bins);
                    break;
                case "mat":
                    readOneFrameMat(currentName, out Az_bins, out fft_bins, out Ang_plot, out RP_dBm);
                    break;
                default:
                    MessageBox.Show("Filetype " + frameType + " not supported yet.");
                    break;
            }

            // detection main algorithm
            bool isFirstFrame = false;
            if (currentId == firstId) isFirstFrame = true;            
            
            if (!justShowOriginal)// && frameType == "mat")
            {
                ranges = new short[16192];
                angles = new short[16192];

                sw.Restart();
                pixNum = InterfaceUWCW.MainInterface(RP_dBm, (short)fft_bins, (short)Az_bins, Ang_plot, Ang_endpt, fft_dR, MagMax, MagMin, isFirstFrame, ranges, angles);
                sw.Stop();
                totalTime += (double)sw.ElapsedMilliseconds / 1000.0;
                totalFrames++;
                showTime();
                
                InterfaceUWCW.GetScoreRange(ref this.score, ref this.prevScore, ref this.minRange, ref this.maxRange);
                file.WriteLine(score);
            }
            
            // paint frame
            Image img = paintFrame(RP_dBm, (int)Az_bins, (int)fft_bins, MagMax, MagMin, ranges, angles, (int)pixNum);
            if (currentId % 2 == 0)
                img.RotateFlip(RotateFlipType.RotateNoneFlipXY);
            else
                img.RotateFlip(RotateFlipType.RotateNoneFlipY);

            // display frame and other parameters, show result, play sound
            dspOriginal.Image = img;
            displayAngleParameters(Az_cnt, Ang_plot);
            playSoundShowLED();
            currentId++;
            cablesShown = true;
        }

        private float[] cnt2Ang(short[] Az_cnt)
        {
            float[] ang = new float[Az_cnt.Length];
            short midCnt = Az_cnt[0];
            midCnt += Az_cnt[Az_cnt.Length - 1];
            midCnt /= 2;

            for (int i = 0; i < Az_cnt.Length; i++)
            {
                ang[i] = (float)(Az_cnt[i] - midCnt) * 0.04f;
            }

            return ang;
        }

        private void playSoundShowLED()
        {
            if (score >= 0.5)
            {
                SystemSounds.Beep.Play();
                for (int i = 0; i < LEDXs.Length; i++)
                {
                    LEDColors[i] = Color.Red;
                }
                txtRange.Text = string.Format("{0:0.0}m ~ {1:0.0}m", minRange, maxRange);
            }
            else if (score >= 0.25 || prevScore >= 0.25)
            {
                LEDColors[0] = Color.LightGreen;
                LEDColors[1] = Color.Orange;
                txtRange.Text = "";
            }
            else
            {
                for (int i = 0; i < LEDXs.Length; i++)
                {
                    LEDColors[i] = Color.LightGreen;
                }
                txtRange.Text = "";
            }

            txtScore.Text = string.Format("{0:0.0}%", score * 100);
            if (currentId == firstId)
                txtPrevScore.Text = string.Format("{0:0.0}%", score * 100);
            else
                txtPrevScore.Text = string.Format("{0:0.0}%", prevScore * 100);

            groupBox1.Refresh();
        }
                
        // read one frame data of the .mat file type
        private void readOneFrameMat(string currentName, out short Az_bins, out int range_bins, out float[] Ang_plot, out float[] RP_dBm)
        {
            Az_bins = 0;
            range_bins = 0;
            Ang_plot = null;
            RP_dBm = null;

            MatFileReader mfr = new MatFileReader(currentName);
            foreach (MLDouble mld in mfr.Data)
            {
                string name = System.Text.ASCIIEncoding.ASCII.GetString(mld.GetNameToByteArray());
                switch (name)
                {
                    case "Ang_plot":
                        Az_bins = (short)mld.Dimensions[1];
                        Ang_plot = new float[Az_bins];
                        for (int i = 0; i < Az_bins; i++)
                        {
                            Ang_plot[i] = (float)mld.GetReal(i);
                        }                        
                        break;
                    case "RP_dBm":
                        range_bins = (short)mld.Dimensions[0];
                        Az_bins = (short)mld.Dimensions[1];
                        RP_dBm = new float[range_bins * Az_bins];
                        int k = 0;
                        for (int i = 0; i < range_bins; i++)
                        {
                            for (int j = 0; j < Az_bins; j++)
                            {
                                RP_dBm[k++] = (float)mld.GetReal(i + j * range_bins);
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }

        // read one frame data of the .bfr file type
        private void readOneFrameBfr(string frameName, out Int16 Az_bins, out Int16[] Az_cnt, out Single[] RP_dBm)
        {
            FileStream fs = null;
            BinaryReader r = null;
            fs = new FileStream(frameName, FileMode.Open, FileAccess.Read);
            r = new BinaryReader(fs);

            // read Az_bins and Az_cnt
            Az_bins = r.ReadInt16();
            Az_cnt = new Int16[Az_bins];
            for (int i = 0; i < Az_bins; i++)
            {
                Az_cnt[i] = r.ReadInt16();
            }

            // read RP_dBm data
            RP_dBm = new Single[Az_bins * fft_bins];
            for (int i = 0; i < Az_bins * fft_bins; i++)
            {
                RP_dBm[i] = r.ReadSingle();
            }

            r.Close();
            fs.Close();
        }     

        private Image paintFrame(float[] RP_dBm, int width, int height, float MagMax, float MagMin, short[] ranges, short[] angles, int pixNum, bool showCables = true)
        {
            Bitmap bmp = new Bitmap(width, height, System.Drawing.Imaging.PixelFormat.Format24bppRgb);
                        
            // lock bits
            Rectangle rect = new Rectangle(0, 0, bmp.Width, bmp.Height);
            System.Drawing.Imaging.BitmapData bmpData = bmp.LockBits(rect, System.Drawing.Imaging.ImageLockMode.WriteOnly, bmp.PixelFormat);

            // get the address of the first pixel
            IntPtr ptr = bmpData.Scan0;

            // declare an arry to hold the bytes of the bitmap
            int bytes = Math.Abs(bmpData.Stride) * bmpData.Height;
            byte[] bgrValues = new byte[bytes]; 

            // set the bgr array values
            colorMapData(RP_dBm, width, height, MagMax, MagMin, bgrValues, Math.Abs(bmpData.Stride));
            int halfLineWidth = 2;
            if (fft_bins > 3000)
            {
                halfLineWidth *= 2;
            }
            if (showCables)
            {
                paintCables(bgrValues, width, height, Math.Abs(bmpData.Stride), ranges, angles, pixNum, halfLineWidth, Color.Red);
            }

            // copy the RGB values to the bitmap
            System.Runtime.InteropServices.Marshal.Copy(bgrValues, 0, ptr, bytes);

            // unlock the bits
            bmp.UnlockBits(bmpData);

            return bmp;
        }

        private void paintCables(byte[] bgrValues, int width, int height, int stride, short[] ranges, short[] angles, int pixNum, int halfLineWidth, Color color)
        {
            byte R = color.R;
            byte G = color.G;
            byte B = color.B;

            int x, y;
            for (int i = 0; i < pixNum; i++)
            {                
                for (int l = -halfLineWidth; l <= halfLineWidth; l++)
                {                   
                    y = (int)(ranges[i]) + l;
                    x = (int)(angles[i]);

                    if (y < 0 || y >= height)
                        continue;
                    if (x < 0 || x >= width)
                        continue;

                    bgrValues[y * stride + 3 * x    ] = B;
                    bgrValues[y * stride + 3 * x + 1] = G;
                    bgrValues[y * stride + 3 * x + 2] = R;
                }
            }
        }

        private void colorMapData(float[] RP_dBm, int width, int height, float MagMax, float MagMin, byte[] bgrValues, int stride)
        {
            int l = 0, k = 0;
            float dataF;
            int dataI;
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    dataF = RP_dBm[l++];
                    dataF -= MagMin;
                    dataF /= (MagMax - MagMin);
                    dataF *= 63;
                    dataI = (int)dataF;
                    if (dataI < 0) dataI = 0;
                    if (dataI > 63) dataI = 63;
                    
                    bgrValues[k + 3 * j  ] = COLORMAP_JET[dataI, 2];
                    bgrValues[k + 3 * j+1] = COLORMAP_JET[dataI, 1];
                    bgrValues[k + 3 * j+2] = COLORMAP_JET[dataI, 0];
                }

                k += stride;
            }
        }

        // flip from column-wise to row-wise
        private float[] flipColToRow(float[] RP_dBm, int width, int height)
        {
            float[] flp = new float[width * height];
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    flp[i * width + j] = RP_dBm[j * height + i];
                }
            }

            return flp;
        }

        // LED paramters and colors
        int LEDSize = 40;
        int[] LEDXs = { 50, 50};
        int[] LEDYs = { 50, 120};
        Color[] LEDColors = { Color.LightGreen, Color.LightGreen, Color.LightGreen };

        private void groupBox1_Paint(object sender, PaintEventArgs e)
        {
            e.Graphics.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
            for (int i = 0; i < LEDXs.Length; i++)
            {
                e.Graphics.DrawEllipse(Pens.Black, LEDXs[i], LEDYs[i], LEDSize, LEDSize);
                e.Graphics.FillEllipse(new SolidBrush(LEDColors[i]), LEDXs[i], LEDYs[i], LEDSize, LEDSize);
            }
        }        

        private void Form1_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (file != null) file.Close();
        }

        private void dspOriginal_Click(object sender, EventArgs e)
        {
            if (!loaded)
            {
                return;
            }

            // paint frame
            Image img = paintFrame(RP_dBm, (int)Az_bins, (int)fft_bins, MagMax, MagMin, ranges, angles, (int)pixNum, !cablesShown);
            cablesShown = !cablesShown;
            if ((currentId-1) % 2 == 0)
                img.RotateFlip(RotateFlipType.RotateNoneFlipXY);
            else
                img.RotateFlip(RotateFlipType.RotateNoneFlipY);

            // display frame and other parameters, show result, play sound
            dspOriginal.Image = img;
        }


        private void justShowOriginalToolStripMenuItem_Click(object sender, EventArgs e)
        {
            justShowOriginal = justShowOriginalToolStripMenuItem.Checked;

            if (justShowOriginal)
            {
                resetVariables();
            }
        }

        private void useOriginalAlgorithmWoTrackingToolStripMenuItem_Click(object sender, EventArgs e)
        {
            useOriginalAlgrm = useOriginalAlgorithmWoTrackingToolStripMenuItem.Checked;

            InterfaceUWCW.SetAlgorithm(useOriginalAlgrm);
        }

        private void showProcessingTimeToolStripMenuItem_Click(object sender, EventArgs e)
        {
            showTime();
        }

        private void showTime()
        {
            if (showProcessingTimeToolStripMenuItem.Checked)
            {
                string s = string.Format("Total time:\n{0:0.000}\nTotal frames:\n{1}\nAverage time:\n{2:0.000}",
                    totalTime, totalFrames, totalTime / (double)totalFrames);
                labelTime.Text = s;
            }
            else
            {
                labelTime.Text = "";
            }
        }
    }
}
