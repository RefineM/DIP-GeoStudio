using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OSGeo.GDAL;
using System.IO;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace dip_1
{
   public static class cls_fusion
    {

        #region brovey变换
        //BROVEY变换
        public static void Brovey(Dataset msidata,Dataset pandata,string savepath,out double[]pingjia)
        {
            cls_basicfunc bc = new cls_basicfunc();
            cls_saveFiles sf = new cls_saveFiles();

            int xsize = pandata.RasterXSize;
            int ysize = pandata.RasterYSize;

            //重采样 以全色为基准
            string path = "C:\\Users\\FAN\\Desktop\\000.tif";
            sf.ChangeFbl(msidata, path,xsize, ysize);
            Dataset msi_cyy = Gdal.Open(path, Access.GA_ReadOnly);
           
            List<int[]> pan = bc.getvalue(pandata);//全色
            List<int[]> msi = bc.getvalue(msi_cyy);//多光谱

            List<int[]> out_msi = new List<int[]>();//多光谱融合结果
            foreach (int[] items in msi)
            { out_msi.Add(new int[xsize * ysize]); }
            
            //多光谱各波段DN均值
            int[]avg_msi =new int[xsize *ysize];
            for (int i=0;i<xsize*ysize;i++)
            {
                for(int band=0;band<msi.Count;band++)
                {
                    avg_msi[i] += msi[band][i];
                }
                avg_msi[i] = avg_msi[i]/msi.Count;
            }

            for (int i = 0; i < xsize * ysize; i++)
            {
                for (int band = 0; band < msi.Count; band++)
                {
                    out_msi[band][i] = msi[band][i]*pan[0][i]/avg_msi[i];
                }
            }
           
            sf.SaveFromDataset(msi_cyy,savepath,out_msi,false);
            bc.deleteFile(path);

            //评价指标
            pingjia = new double[6];

            //RGB灰度化
            cls_segment sg = new cls_segment();
            int[] bef_gray = sg.ToGray(msi[0], msi[1], msi[2]);
            int[] aft_gray = sg.ToGray(out_msi[0], out_msi[1], out_msi[2]);
            
            pingjia[0] = Entropy(bef_gray);  //融合前多光谱图像信息熵
            pingjia[1] = Entropy(aft_gray);  //融合后多光谱图像信息熵

            pingjia[2] = coEntropy(bef_gray,aft_gray); //融合前后多光谱图像联合熵

            pingjia[3] = defination(bef_gray,xsize,ysize); //融合前多光谱图像清晰度（平均梯度）
            pingjia[4] = defination(aft_gray,xsize,ysize); //融合后多光谱图像清晰度

            pingjia[5] = fidelity(bef_gray,aft_gray);  //融合后多光谱图像逼真度（偏差）

        }
        #endregion
       
        #region IHS变换
        //IHS变换
        public static void IHS(Dataset msidata, Dataset pandata, string savepath,out double[] pingjia)
        {
            cls_basicfunc bc = new cls_basicfunc();
            cls_saveFiles sf = new cls_saveFiles();

            int xsize = pandata.RasterXSize;
            int ysize = pandata.RasterYSize;

            //以全色为基准重采样msi
            string path = "C:\\Users\\FAN\\Desktop\\000.tif";
            sf.ChangeFbl(msidata, path, xsize, ysize);
            Dataset msi_cyy = Gdal.Open(path, Access.GA_ReadOnly);

            //全色 多光谱 rgb
            List<int[]> pan = bc.getvalue(pandata);
            double[] pan_I = new double[pan[0].Length];
            for(int i=0;i< pan[0].Length; i++)
            {
                pan_I[i] = pan[0][i] * 1.0;
            }
            List<int[]> msi = bc.getvalue(msi_cyy);
            List<int[]> rgb = new List<int[]>();
            for (int i=0;i<3;i++)
            {
                rgb.Add(msi[i]);
            }

            //多光谱融合结果
            List<int[]> outmsi = new List<int[]>();
            foreach (int []items in rgb)
            { outmsi.Add(new int[xsize * ysize]); }

            //step1:把rgb空间转换为ihs空间
            List<double[]> msi_IHS = RGBtoIHS(rgb);
           
            //step2:将pan_I和msi_I直方图匹配
            double[] msi_I = msi_IHS[0];
            double[] res_I = histoMatch(pan_I,msi_I);
            List<double[]> msi_I2HS = new List<double[]>() { res_I , msi_IHS[1] , msi_IHS[2] };

            //step4:把ihs空间转换为rgb空间
            List<int[]> msi_rgb = IHStoRGB(msi_I2HS);

            //保存
            sf.SaveFromDataset(msi_cyy, savepath, msi_rgb, false);

            //评价指标
            pingjia = new double[6];

            //RGB灰度化
            cls_segment sg = new cls_segment();
            int[] bef_gray = sg.ToGray(msi[0], msi[1], msi[2]);
            int[] aft_gray = sg.ToGray(msi_rgb[0], msi_rgb[1], msi_rgb[2]);

            pingjia[0] = Entropy(bef_gray);  //融合前多光谱图像信息熵
            pingjia[1] = Entropy(aft_gray);  //融合后多光谱图像信息熵

            pingjia[2] = coEntropy(bef_gray, aft_gray); //融合前后多光谱图像联合熵

            pingjia[3] = defination(bef_gray, xsize, ysize); //融合前多光谱图像清晰度（平均梯度）
            pingjia[4] = defination(aft_gray, xsize, ysize); //融合后多光谱图像清晰度

            pingjia[5] = fidelity(bef_gray, aft_gray);  //融合后多光谱图像逼真度（偏差）
        }
        // 使用“几何推导法”进行颜色空间转换
        // RGB -> IHS
        public static List<double[]> RGBtoIHS (List<int[]>RGB)
        {
            int[] r = RGB[0];
            int[] g = RGB[1];
            int[] b = RGB[2];
       
           // 先对RGB归一化到[0，1]
            int rmin = minvalue(r);
            int gmin = minvalue(g);
            int bmin = minvalue(b);
            int rmax = maxvalue(r);
            int gmax = maxvalue(g);
            int bmax = maxvalue(b);
            double[] normR = new double[r.Length];
            double[] normG = new double[r.Length];
            double[] normB = new double[r.Length];
            for(int i=0;i<r.Length;i++)
            { 
                normR[i] = (r[i]-rmin)*1.0 / (rmax - rmin);
                normG[i] = (g[i]-gmin)*1.0 / (gmax - gmin);
                normB[i] = (b[i]-bmin)*1.0 / (bmax - bmin);
            }

            // 注意IHS的数值范围: (1)I和S本身在[0,1]; (2) H在[0,2*pi],要归一化到[0，1]
            double[] I = new double[r.Length];
            double[] H = new double[r.Length];
            double[] S = new double[r.Length];
            double[]sita = new double[r.Length];
            for (int i=0;i<r.Length;i++)
            {
                double[] rgb = { normR[i], normG[i], normB[i]};
                I[i] = (normR[i] + normG[i] + normB[i]) /3;
                S[i] = 1 - 3 * minvalueD(rgb) / (normR[i] + normG[i] + normB[i]);
                sita[i] = Math.Acos(0.5*(2* normR[i]- normG[i]- normB[i])/
                    Math.Sqrt(Math.Pow((normR[i]- normG[i]),2)+(normR[i]- normB[i])*(normG[i]- normB[i])));
                if (normG[i] >= normB[i])
                {
                    H[i] = sita[i]/ (2 * Math.PI);
                }else
                {
                    H[i] = (2*Math.PI-sita[i])/ (2 * Math.PI);
                }
            }
            List<double[]> IHS = new List<double[]>() {I,H,S};
            return IHS;
        }
        // RGB <- IHS
        public static List<int[]> IHStoRGB(List<double[]> IHS)
        {
            // 先将H范围扩充到[0-2*pi]
            double[] H = new double[IHS[1].Length];
            for(int i=0;i<IHS[1].Length;i++)
            {
                H[i] = IHS[1][i] * 2 * Math.PI;
            }
            double[] I = IHS[0];
            double[] S = IHS[2];

            int[] R = new int[I.Length];
            int[] G = new int[I.Length];
            int[] B = new int[I.Length];

            for (int i = 0; i < I.Length; i++)
            {
                if(H[i]>=0 && H[i]<(120*Math.PI/180))
                {
                    R[i] = (int)(I[i] * (1 + S[i] * Math.Cos(H[i]) / Math.Cos(60 * Math.PI / 180 - H[i])));
                    B[i] = (int)(I[i] * (1 - S[i]));
                    G[i] = (int)(3* I[i] - R[i] - B[i]);
                }
                if(H[i] >= (120 * Math.PI / 180) && H[i] < (240 * Math.PI / 180))
                {
                    G[i] = (int)(I[i] * (1 + S[i] * Math.Cos(H[i]-120*Math.PI/180) / Math.Cos(180 * Math.PI / 180 - H[i])));
                    R[i] = (int)(I[i] * (1 - S[i]));
                    B[i] = (int)(3 * I[i] - R[i] - G[i]);
                }
                if(H[i] >= (240 * Math.PI / 180) && H[i] < (360 * Math.PI / 180))
               {
                    B[i] = (int)(I[i] * (1 + S[i] * Math.Cos(H[i] -240 * Math.PI / 180) / Math.Cos(300 * Math.PI / 180 - H[i])));
                    G[i] = (int)(I[i] * (1 - S[i]));
                    R[i] = (int)(3 * I[i] - G[i] - B[i]);
                }
            }
            List<int[]> RGB = new List<int[]>() {R,G,B};
            return RGB;
        }
        // minvalue
        public static int minvalue(int[]rgb)
        {
            int min = rgb[0];
            for(int i=1;i<rgb.Length;i++)
            {
                if (rgb[i] < min)
                    min = rgb[i];
            }
            return min;
        }
        // maxvalue
        public static int maxvalue(int[] rgb)
        {
            int max = rgb[0];
            for (int i = 1; i < rgb.Length; i++)
            {
                if (rgb[i] > max)
                    max = rgb[i];
            }
            return max;
        }
        public static double minvalueD(double[] rgb)
        {
            double min = rgb[0];
            for (int i = 1; i < rgb.Length; i++)
            {
                if (rgb[i] < min)
                    min = rgb[i];
            }
            return min;
        }
        // maxvalue
        public static double maxvalueD(double[] rgb)
        {
            double max = rgb[0];
            for (int i = 1; i < rgb.Length; i++)
            {
                if (rgb[i] > max)
                    max = rgb[i];
            }
            return max;
        }

        //直方图匹配(两个波段间)
        public static double[] histoMatch(double[] pan,double[] msi)
        {
            //匹配后的像素值集合
            double[] res = new double[msi.Length];
            double min = minvalueD(msi);//归一化前msi的最值
            double max = maxvalueD(msi);

            //直方图均衡化
            double[] balpan = histoBalance(pan, out int[]normpan, out int[] mappan);
            double[] balmsi = histoBalance(msi, out int[]normmsi, out int[] mapmsi);

            //求均衡化之后的累计概率分布直方图
            double[] miduhistoPan = leijimiduhisto(balpan);
            double[] miduhistoMsi = leijimiduhisto(balmsi);

            // 求出msi累积直方图每一个灰度级的概率和pan累积直方图每一个灰度级概率之差
            double[,] diff = new double [miduhistoMsi.Length,miduhistoPan.Length];
            for(int i=0;i<miduhistoMsi.Length;i++)
            {
                for(int j=0;j<miduhistoPan.Length;j++)
                {
                    diff[i, j] = Math.Abs(miduhistoMsi[i] - miduhistoPan[j]);
                }
            }

            // 求出概率差最小时msi对应pan的灰度级
            // id[i]=j i是balmsi的灰度级 j是balpan的灰度级
            int[] mapbal = new int[miduhistoMsi.Length];
            double minv = 0.0;
            for (int i = 0; i < miduhistoMsi.Length; i++)
            {
                minv = diff[i, 0];
                mapbal[i] = 0;
                for (int j = 1; j < miduhistoPan.Length; j++)
                {
                    if (diff[i, j] < minv)
                    {
                        minv = diff[i, j];
                        mapbal[i] = j;
                    }
                } 
            }

            // 映射
            for (int i = 0; i < msi.Length; i++)
            {
                  res[i] = mapbal[mapmsi[normmsi[i]]];
            }
            return res;
        }

        //直方图均衡化(单波段，灰度级未归算)
        public static double[] histoBalance(double[] orivalue,out int[] norm,out int[]map )
        {
            //pre:准备
            int[] histo = new int[256];//频数直方图
            double[] miduhisto = new double[256];//累计概率直方图
            map = new int[256];//查找表
      
            double min = minvalueD(orivalue);//归一化前的最值
            double max = maxvalueD(orivalue);
         
            norm = new int[orivalue.Length];//归一化后像元值
            double[] res = new double[orivalue.Length];//处理后的像素值集合
        
            //step1：【归一化】 把像元值归到0-255区间
            for (int i = 0; i < orivalue.Length; i++)
            {
                norm[i] = (int)((orivalue[i] - min) * 255 / (max - min));
            }
         
            //step2:【直方图均衡化】
            // 频数直方图
            for (int i = 0; i < orivalue.Length; i++) histo[norm[i]]++;
            // 累计概率直方图
            for (int i = 0; i < 256; i++)
            {
                if (i == 0)
                {
                    miduhisto[i] = histo[0];
                }
                else
                {
                    miduhisto[i] = miduhisto[i - 1] + histo[i];
                }
            }
            for (int i = 0; i < 256; i++)
            {
                miduhisto[i] /= orivalue.Length;
            }
            // 均衡化操作
            for (int i = 0; i < 256; i++)
            {
                map[i] = (int)(255 * miduhisto[i] + 0.5);
            }
            // 映射
            for (int i = 0; i < orivalue.Length; i++)
            {
                res[i] = map[norm[i]];
            }
            return res;
        }

        //求累计概率密度分布直方图（单波段，灰度级已归算到0-255）
        public static double[] leijimiduhisto(double[] orivalue)
        {
            double[] histo = new double[256];
            double[] miduhisto = new double[256];
            // 频数直方图
            for (int i = 0; i < orivalue.Length; i++)
            {
                histo[(int)(orivalue[i])]++;
            }
            // 累计概率直方图
            for (int i = 0; i < 256; i++)
            {
                if (i == 0)
                {
                    miduhisto[i] = histo[0];
                }
                else
                {
                    miduhisto[i] = miduhisto[i - 1] + histo[i];
                }
            }
            for (int i = 0; i < 256; i++)
            {
                miduhisto[i] /= orivalue.Length;
            }
            return miduhisto;
        }
        #endregion

        #region 融合评价
        //信息熵
        public static double Entropy(int[] band)
        {
            //归算到0-255
            cls_segment sg = new cls_segment();
            int[] norband = sg.ToNormal(band);

            //求概率分布直方图
            double[] histo = sg.GetNormalHisto(norband);
            double[] value = new double[256];
            double entropy = 0.0;

            //求信息熵
            for (int i = 0; i < 256; i++)
            {
                if(histo[i]!=0)
                {
                    value[i] = histo[i] * Math.Log(1.0 / histo[i], 2.0);
                }else
                {
                    value[i] = 0.0;
                }
            }

            for(int i=0;i<256;i++)
            {
                entropy += value[i];
            }
            return entropy;
        }

        //联合熵
        public static double coEntropy(int[] band1, int[] band2)
        {
            //归0-255
            cls_segment sg = new cls_segment();
            int[] norband1 = sg.ToNormal(band1);
            int[] norband2 = sg.ToNormal(band2);

            //求联合概率分布直方图
            int num = norband1.Length;
            double[,] cohisto = new double[256, 256];
            double entropy = 0.0;
            for (int i=0;i<num;i++)
            {
               cohisto[norband1[i],norband2[i]]+=1.0;        
            }
            for (int i = 0; i < 256; i++)
            {
                for (int j = 0; j < 256; j++)
                {
                    cohisto[i,j]/= num;
                }
            }
            //求联合熵
            for (int i = 0; i < 256; i++)
            {
                for (int j = 0; j < 256; j++)
                {
                    if(cohisto[i,j]!=0)
                    {
                        entropy += cohisto[i, j] * Math.Log(1/cohisto[i, j], 2);
                    }
                }
            }
            return entropy;
        }

        //清晰度 (平均梯度)
        public static double defination(int[] band, int xsize, int ysize)
        {
            double defination = 0.0;
            for (int i = 1; i < xsize - 1; i++)
                for (int j = 1; j < ysize - 1; j++)
                {
                    defination += Math.Sqrt(
                                 ( Math.Pow((band[i + j * xsize] - band[i + (j + 1) * xsize]), 2) +
                                   Math.Pow((band[i + j * xsize] - band[(i + 1) + j * xsize]), 2)
                                  ) * 0.5);
                }
            defination = defination / (xsize - 1) / (ysize - 1);
            return defination;
        }

        //逼真度（偏差系数）
        public static double fidelity(int[] band1, int[] band2)
        {
            double f = 0.0;
            for (int i = 0; i < band1.Length; i++)
            {
                f += Math.Abs(band1[i] - band2[i]) / band1[i];
            }
            f = f / band1.Length;
            return f;
        }
        #endregion

        #region IHS(2)
        private static List<int[]> GetDNs(Dataset ds, int xSize, int ySize)
        {
            List<int[]> dns = new List<int[]>();
            for (int i = 1; i <= ds.RasterCount; i++)
            {
                int[] temp = new int[xSize * ySize];
                double[] mm = { 0, 0 };
                ds.GetRasterBand(i).ComputeRasterMinMax(mm, 0);
                ds.GetRasterBand(i).ReadRaster(0, 0, ds.RasterXSize, ds.RasterYSize, temp, xSize, ySize, 0, 0);
                for (int j = 0; j < xSize * ySize; j++)
                {
                    temp[j] = (int)(((temp[j] - mm[0]) / (mm[1] - mm[0])) * 255);
                }
                dns.Add(temp);
            }
            return dns;
        }
        /// <summary>
        /// IHS变换
        /// </summary>
        /// <param name="msi">多光谱影像数据</param>
        /// <param name="pan">全色影像数据</param>
        /// <param name="filePath">保存路径</param>
        public static void IHS2(Dataset msi, Dataset pan, string filePath, out double[] pingjia)
        {
            // 全色影像单波段数据
            int xSize = pan.RasterXSize;
            int ySize = pan.RasterYSize;
            int[] pan_r = GetDNs(pan, xSize, ySize)[0];

            //以全色为基准重采样msi
            string path = "C:\\Users\\FAN\\Desktop\\000.tif";
            cls_saveFiles sf = new cls_saveFiles();
            sf.ChangeFbl(msi, path, xSize, ySize);
            Dataset msi_cyy = Gdal.Open(path, Access.GA_ReadOnly);

            // 多光谱影像多波段数据
            List<int[]> msi_all = GetDNs(msi_cyy, xSize, ySize);
            int[] msi_r = msi_all[0];
            int[] msi_g = msi_all[1];
            int[] msi_b = msi_all[2];

            // RGB空间变换到IHS空间
            List<double[]> ihs = RecRGBToIHS(new List<int[]>() { msi_r, msi_g, msi_b });

            // 对全色影像和IHS空间中的亮度分量I进行直方图匹配
            int[] matched = HistoMatch(pan_r, ihs[0].ToInt());

            // 用全色影像I'代替IHS空间的亮度分量
            ihs[0] = matched.ToDouble();

            // 将I'HS逆变换到RGB空间
            List<int[]> band = RecIHSToRGB(ihs);
            cls_saveFiles sa = new cls_saveFiles();
            sa.SaveFromDataset(pan, filePath, band,false);

            //评价指标
            pingjia = new double[6];

            //RGB灰度化
            cls_segment sg = new cls_segment();
            int[] bef_gray = sg.ToGray(msi_all[0],msi_all[1],msi_all[2]);
            int[] aft_gray = sg.ToGray(band[0], band[1], band[2]);

            pingjia[0] = Entropy(bef_gray);  //融合前多光谱图像信息熵
            pingjia[1] = Entropy(aft_gray);  //融合后多光谱图像信息熵

            pingjia[2] = coEntropy(bef_gray, aft_gray); //融合前后多光谱图像联合熵

            pingjia[3] = defination(bef_gray, xSize, ySize); //融合前多光谱图像清晰度（平均梯度）
            pingjia[4] = defination(aft_gray, xSize, ySize); //融合后多光谱图像清晰度

            pingjia[5] = fidelity(bef_gray, aft_gray);  //融合后多光谱图像逼真度（偏差）
        }

        private static List<double[]> RecRGBToIHS(List<int[]> rgb)
        {
            int cou = rgb[0].Count();
            double[] i = new double[cou];
            double[] h = new double[cou];
            double[] s = new double[cou];

            for (int j = 0; j < cou; j++)
            {
                double r = rgb[0][j] * 1.0;
                double g = rgb[1][j] * 1.0;
                double b = rgb[2][j] * 1.0;
                double _i = r + g + b;

                i[j] = _i / 3.0;

                double min = Math.Min(Math.Min(r, g), b) * 1.0;
                if (r == g && g == b)
                {
                    h[j] = 0;
                    s[j] = 0;
                }
                else
                {
                    if (min == b)
                    {
                        h[j] = (g - b) / (_i - 3.0 * b);
                        s[j] = (_i - 3.0 * b) / _i;
                    }
                    else if (min == r)
                    {
                        h[j] = (b - r) / (_i - 3.0 * r) + 1;
                        s[j] = (_i - 3.0 * r) / _i;
                    }
                    else
                    {
                        h[j] = (r - g) / (_i - 3.0 * g) + 2;
                        s[j] = (_i - 3.0 * g) / _i;
                    }
                }
            }


            return new List<double[]>() { i, h, s };
        }

        private static List<int[]> RecIHSToRGB(List<double[]> ihs)
        {
            int cou = ihs[0].Count();
            int[] r = new int[cou];
            int[] g = new int[cou];
            int[] b = new int[cou];

            for (int j = 0; j < cou; j++)
            {
                double i = ihs[0][j] * 1.0;
                double h = ihs[1][j] * 1.0;
                double s = ihs[2][j] * 1.0;

                if (h >= 2.0 && h < 3.0)
                {
                    r[j] = (int)(i * (1.0 - 7 * s - 3 * s * h));
                    g[j] = (int)(i * (1.0 - s));
                    b[j] = (int)(i * (1.0 + 8 * s - 2 * s * h));
                }
                else if (h >= 1.0 && h < 2.0)
                {
                    r[j] = (int)(i * (1.0 - s));
                    g[j] = (int)(i * (1.0 + 5 * s - 3 * s * h));
                    b[j] = (int)(i * (1.0 - 4 * s + 3 * s * h));
                }
                else
                {
                    r[j] = (int)(i * (1.0 + 2 * s - 3 * s * h));
                    g[j] = (int)(i * (1.0 - 1 * s + 3 * s * h));
                    b[j] = (int)(i * (1.0 - s));
                }

                if (r[j] > 255) r[j] = 255;
                if (g[j] > 255) g[j] = 255;
                if (b[j] > 255) b[j] = 255;
                if (r[j] < 0) r[j] = 0;
                if (g[j] < 0) g[j] = 0;
                if (b[j] < 0) b[j] = 0;
            }

            return new List<int[]>() { r, g, b };
        }

        /// <summary>
        /// int型数组转为double型
        /// </summary>
        /// <param name="arr"></param>
        /// <returns></returns>
        public static int[] ToInt(this double[] arr)
        {
            int[] fin = new int[arr.Length];
            for (int i = 0; i < arr.Length; i++)
            {
                fin[i] = (int)arr[i];
            }
            return fin;
        }

        /// <summary>
        /// double型数组转为int型
        /// </summary>
        /// <param name="arr"></param>
        /// <returns></returns>
        public static double[] ToDouble(this int[] arr)
        {
            double[] fin = new double[arr.Length];
            for (int i = 0; i < arr.Length; i++)
            {
                fin[i] = arr[i];
            }
            return fin;
        }

        /// <summary>
        /// 直方图匹配
        /// </summary>
        /// <param name="pan">多光谱影像</param>
        /// <param name="res">高分辨率影像</param>
        /// <returns></returns>
        public static int[] HistoMatch(int[] pan, int[] res)
        {
            int[] _pan = GetHistograph(pan);
            int[] _res = GetHistograph(res);

            double[] pan_dn = new double[buckets];
            double[] res_dn = new double[buckets];
            for (int i = 0; i < buckets; i++)
            {
                pan_dn[i] = ((_pan[i] * 1.0) / pan.Length);
                res_dn[i] = ((_res[i] * 1.0) / res.Length);
            }
            for (int i = 1; i < buckets; i++)
            {
                pan_dn[i] = pan_dn[i] + pan_dn[i - 1];
                res_dn[i] = res_dn[i] + res_dn[i - 1];
            }
            double diffAR = 0.0, diffBR = 0.0;
            byte kR = 0;
            byte[] mapPixelR = new byte[buckets];
            for (int i = 0; i < buckets; i++)
            {
                diffBR = 1;
                for (int j = kR; j < buckets; j++)
                {
                    diffAR = Math.Abs(pan_dn[i] - res_dn[j]);
                    if (diffAR - diffBR < 1.0E-08)
                    {
                        diffBR = diffAR;
                        kR = (byte)j;
                    }
                    else
                    {
                        kR = (byte)Math.Abs(j - 1);
                        break;
                    }
                }
                if (kR == buckets - 1)
                {
                    for (int l = i; l < buckets; l++)
                    {
                        mapPixelR[l] = kR;
                    }
                    break;
                }
                mapPixelR[i] = kR;
            }

            int[] fin = new int[res.Length];

            for (int i = 0; i < res.Length; i++)
            {
                fin[i] = mapPixelR[pan[i]];
            }
            return fin;
        }

        /// <summary>
        /// 颜色宽度
        /// </summary>
        public static int buckets = 4096;

        /// <summary>
        /// 获取int型灰度直方图
        /// </summary>
        /// <param name="arr"></param>
        /// <returns></returns>
        public static int[] GetHistograph(int[] arr)
        {
            int[] histo = new int[buckets];
            for (int i = 0; i < arr.Length; i++)
            {
                histo[arr[i]]++;
            }
            return histo;
        }
     
        #endregion

    }
}
