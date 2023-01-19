using System;
using System.Collections.Generic;
using System.Linq;
using OSGeo.GDAL;

namespace dip_1
{
    class cls_segment
    {
        // RGB转Gray
        public int[] ToGray(int[] R, int[] G, int[] B)
        {
            int[] gray = new int[R.Length];
            for (int i = 0; i < gray.Length; i++)
            {
                gray[i] = (int)(0.299 * R[i] + 0.587 * G[i] + 0.114 * B[i]);
            }
            return gray;
        }

        // 归化到[0，255]区间
        public int[] ToNormal(int[] array)
        {
            int min = array.Min();
            int max = array.Max();
            int[] res = new int[array.Length];
            for (int i = 0; i < array.Length; i++)
            {
                res[i] = (array[i] - min) * 255 / (max - min);
            }
            return res;
        }

        // 概率密度直方图
        public double[] GetNormalHisto(int[] array)
        {
            double[] histoGray = new double[256];
            double[] glhistoGray = new double[256];
            int pixels = array.Length;


            for (int i = 0; i < pixels; i++)
            {
                histoGray[array[i]]++;
            }

            for (int i = 0; i < 256; i++)
            {
                glhistoGray[i] = histoGray[i] / pixels;
            }

            return glhistoGray;
        }

        // 双阈值大津法
        public void OtsuDouble(Dataset dt, string path, out int T1, out int T2)
        {
            cls_basicfunc bc = new cls_basicfunc();
            List<int[]> data = bc.getvalue(dt);

            int xsize = dt.RasterXSize;
            int ysize = dt.RasterYSize;
            int nums = xsize * ysize;
            int [] gray = new int[nums];

            if(dt.RasterCount>=3)
            {
                int[] R = data[0];
                int[] G = data[1];
                int[] B = data[2];
                gray = ToGray(R, G, B);//灰度化
                gray = ToNormal(gray);//归一化
            }
            if(dt.RasterCount==1)
            {
                gray = data[0];//灰度化
                gray = ToNormal(gray);//归一化
            }

            double[] histoGray = GetNormalHisto(gray);//概率密度直方图

            T1 = 0;
            T2 = 0;
            double var_max = 0.0;
            for (int t1 = 1; t1 <= 253; t1++)
            {
                double p1 = 0.0, avg1 = 0.0;
                //统计 0<=DN<t1 范围的像元平均值、占比
                for (int i = 0; i < t1; i++)
                {
                    p1 += histoGray[i];
                    avg1 += histoGray[i] * i;
                }
                avg1 = avg1 / p1;
                for (int t2 = t1 + 1; t2 <= 255; t2++)
                {
                    double p2 = 0.0, p3 = 0.0, var = 0.0, avg2 = 0.0, avg3 = 0.0;
                    //统计 t1<= DN <t2 范围的像元平均值、占比
                    for (int i = t1; i < t2; i++)
                    {
                        avg2 += histoGray[i] * i;
                        p2 += histoGray[i];
                    }
                    avg2 = avg2 / p2;
                    //统计 t2 <= DN <=255 范围的像元平均值、占比
                    for (int i = t2; i <= 255; i++)
                    {
                        avg3 += histoGray[i] * i;
                        p3 += histoGray[i];
                    }
                    avg3 = avg3 / p3;
                    //总方差
                    double avg0 = p1 * avg1 + p2 * avg2 + p3 * avg3; 
                    var = p1 * Math.Pow((avg1 - avg0), 2) + p2 * Math.Pow((avg2 - avg0), 2) + p3 * Math.Pow((avg3 - avg0), 2);
                    if (var > var_max)
                    {
                        var_max = var;
                        T1 = t1;
                        T2 = t2;
                    }
                }
            }

            //映射
            for (int i = 0; i < xsize * ysize; i++)
            {
                if (gray[i] < T1)
                {
                    gray[i] = 0;
                }

                if (gray[i] >= T1 && gray[i] <= T2)
                {
                    gray[i] = 127;
                }

                if (gray[i] > T2)
                {
                    gray[i] = 255;
                }
            }

            List<int[]> band = new List<int[]>();
            band.Add(gray);
            cls_saveFiles sa = new cls_saveFiles();
            sa.SaveFromDataset(dt, path, band, false);

        }

        // 单阈值大津法
        public void OtsuSingle(Dataset dt, string path, out int T1)
        {
            cls_basicfunc bc = new cls_basicfunc();
            List<int[]> data = bc.getvalue(dt);

            int xsize = dt.RasterXSize;
            int ysize = dt.RasterYSize;
            int nums = xsize * ysize;
            int[] gray = new int[nums];

            if (dt.RasterCount >= 3)
            {
                int[] R = data[0];
                int[] G = data[1];
                int[] B = data[2];
                gray = ToGray(R, G, B);//灰度化
                gray = ToNormal(gray);//归一化
            }
            if (dt.RasterCount == 1)
            {
                gray = data[0];//灰度化
                gray = ToNormal(gray);//归一化
            }

            double[] histoGray = GetNormalHisto(gray);//概率密度直方图
            T1 = 0;
            double var_max = 0.0;

            for (int t1 = 1; t1 <= 254; t1++)
            {
                double p1 = 0.0, avg1 = 0.0, p2 = 0.0, avg2 = 0.0, var = 0.0;

                for (int i = 0; i < t1; i++)
                {
                    p1 += histoGray[i];
                    avg1 += histoGray[i] * i;
                }
                avg1 = avg1 / p1;

                for (int i = t1; i < 256; i++)
                {
                    p2 += histoGray[i];
                    avg2 += histoGray[i] * i;
                }
                avg2 = avg2 / p2;

                double avg0 = p1 * avg1 + p2 * avg2;
                var = p1 * Math.Pow((avg1 - avg0), 2) + p2 * Math.Pow((avg2 - avg0), 2);
                if (var > var_max)
                {
                    var_max = var;
                    T1 = t1;
                }
            }

            //映射
            for (int i = 0; i < xsize * ysize; i++)
            {
                if (gray[i] < T1)
                {
                    gray[i] = 0;
                }
                else
                {
                    gray[i] = 255;
                }
            }

            List<int[]> band = new List<int[]>();
            band.Add(gray);
            cls_saveFiles sa = new cls_saveFiles();
            sa.SaveFromDataset(dt, path, band, false);
        }

        // 迭代阈值法
        public void iteration(Dataset dt, string path, out int T)
        {
            cls_basicfunc bc = new cls_basicfunc();
            List<int[]> data = bc.getvalue(dt);
            int xsize = dt.RasterXSize;
            int ysize = dt.RasterYSize;
            int nums = xsize * ysize;

            int[] R = data[0]; int[] G = data[1];int[] B = data[2];

            int[] _gray = ToGray(R, G, B);//灰度化
            int[] gray = ToNormal(_gray);//归一化
            double[] histoGray = GetNormalHisto(gray);//概率密度直方图
            double avg0 = gray.Average();//整体平均值
            double t0 = avg0;

            while (1 == 1)
            {
                double p1 = 0.0, avg1 = 0.0, p2 = 0.0, avg2 = 0.0;
                for (int i = 0; i < (int)(t0); i++)
                {
                    p1 += histoGray[i];
                    avg1 += histoGray[i] * i;
                }
                avg1 = avg1 / p1;

                for (int i = (int)(t0); i < 256; i++)
                {
                    p2 += histoGray[i];
                    avg2 += histoGray[i] * i;
                }
                avg2 = avg2 / p2;

                int t1 = (int)((avg1 + avg2) / 2);

                if (Math.Abs(t1 - t0) <= 0.01)
                {
                    T = t1;
                    break;
                }
                else
                {
                    t0 = t1;
                }
            }

            //映射
            for (int i = 0; i < xsize * ysize; i++)
            {
                if (gray[i] < T)
                {
                    gray[i] = 0;
                }
                else
                {
                    gray[i] = 255;
                }
            }

            List<int[]> band = new List<int[]>();
            band.Add(gray);
            cls_saveFiles sa = new cls_saveFiles();
            sa.SaveFromDataset(dt, path, band, false);

        }

        //差值法变化检测
        public void changeDetection(Dataset dt1, Dataset dt2, string savepath, out int T1)
        {
            cls_basicfunc bc = new cls_basicfunc();
            List<int[]> data1 = bc.getvalue(dt1);
            List<int[]> data2 = bc.getvalue(dt2);
            int bands = dt1.RasterCount;
            int pixels = dt1.RasterXSize * dt2.RasterYSize;

            List<int[]> diff = new List<int[]>();
            foreach(int[] i in data1)
            {
                diff.Add(new int[pixels]);
            }

            //求差
            for (int i = 0; i < bands; i++)
            {
                for (int j = 0; j < pixels; j++)
                {
                    diff[i][j] = Math.Abs(data1[i][j] - data2[i][j]);
                }
            }

            cls_saveFiles sa = new cls_saveFiles();
            string p = "C:\\Users\\FAN\\Desktop\\diff.tif";
            sa.SaveFromDataset(dt1, p, diff, false);
            Dataset dif = Gdal.Open(p, Access.GA_ReadOnly);

            //OtsuSingle(dif,savepath, out T1);
            iteration(dif, savepath, out T1);

        }
    }
}

