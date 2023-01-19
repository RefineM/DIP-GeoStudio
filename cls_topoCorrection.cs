using OSGeo.GDAL;
using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;

namespace dip_1
{
    class cls_topoCorrection
    {

        #region 坡度坡向
        //模板填充 返回填充完毕的dem高程数据
        public double[] fulldem(Dataset dt)
        {
            int xSize = dt.RasterXSize;
            int ySize = dt.RasterYSize;
            double[] oldvalue = new double[xSize * ySize];

            Band band = dt.GetRasterBand(1);//dem高程数据为单波段，存储高程信息
            band.ReadRaster(0, 0, xSize, ySize, oldvalue, xSize, ySize, 0, 0);//数据读入到数组中

            //将原图像外扩充一圈，便于坡度计算
            double[] newvalue = new double[(xSize + 2) * (ySize + 2)];
            for (int i = 0; i < xSize; i++)//填充扩充的上下两行
            {
                newvalue[i + 1] = oldvalue[i];
                newvalue[(ySize + 1) * (xSize + 2) + 1 + i] = oldvalue[(ySize - 1) * xSize + i];
            }
            for (int j = 0; j < ySize; j++)//填充扩充的左右两行
            {
                newvalue[(j + 1) * (xSize + 2)] = oldvalue[j * xSize];
                newvalue[(j + 1) * (xSize + 2) + (xSize + 1)] = oldvalue[j * xSize + (xSize - 1)];
            }

            //填充四个角点
            newvalue[0] = oldvalue[0];//左上角点
            newvalue[xSize + 1] = oldvalue[xSize - 1];//右上角点
            newvalue[(xSize + 2) * (ySize + 1)] = oldvalue[xSize * (ySize - 1)];//左下角点；
            newvalue[(xSize + 2) * (ySize + 2) - 1] = oldvalue[xSize * ySize - 1];//右下角点

            //复制原有的像元值
            for (int i = 0; i < dt.RasterXSize; i++)
            {
                for (int j = 0; j < dt.RasterYSize; j++)
                {
                    newvalue[(i + 1) + j * (dt.RasterXSize + 2)] = oldvalue[i + j * dt.RasterXSize];
                }
            }
            return newvalue;
        }
        //模板运算坡度、坡向
        public void OperatorEnforce(Dataset dt, double[,] cons, int[,] ope_X, int[,] ope_Y, double fbl_x,double fbl_y,out int[]slope_res, out int[] aspect_res)
        {
            slope_res = new int[dt.RasterXSize*dt.RasterYSize];
            aspect_res = new int[dt.RasterXSize * dt.RasterYSize];

            for (int i = 1; i < dt.RasterXSize + 1; i++)
            {
                for (int j = 1; j < dt.RasterYSize + 1; j++)
                {
                    double slope_X = 0;
                    slope_X += cons[i - 1, j - 1] * ope_X[0, 0];
                    slope_X += cons[i, j - 1] * ope_X[0, 1];
                    slope_X += cons[i + 1, j - 1] * ope_X[0, 2];
                    slope_X += cons[i - 1, j] * ope_X[1, 0];
                    slope_X += cons[i, j] * ope_X[1, 1];
                    slope_X += cons[i + 1, j] * ope_X[1, 2];
                    slope_X += cons[i - 1, j + 1] * ope_X[2, 0];
                    slope_X += cons[i, j + 1] * ope_X[2, 1];
                    slope_X += cons[i + 1, j + 1] * ope_X[2, 2];
                    double slope_Y = 0;
                    slope_Y += cons[i - 1, j - 1] * ope_Y[0, 0];
                    slope_Y += cons[i, j - 1] * ope_Y[0, 1];
                    slope_Y += cons[i + 1, j - 1] * ope_Y[0, 2];
                    slope_Y += cons[i - 1, j] * ope_Y[1, 0];
                    slope_Y += cons[i, j] * ope_Y[1, 1];
                    slope_Y += cons[i + 1, j] * ope_Y[1, 2];
                    slope_Y += cons[i - 1, j + 1] * ope_Y[2, 0];
                    slope_Y += cons[i, j + 1] * ope_Y[2, 1];
                    slope_Y += cons[i + 1, j + 1] * ope_Y[2, 2];
                    //slopeAndAspect
                    double slope = Math.Sqrt(Math.Pow(slope_X / (8 * fbl_x), 2) + Math.Pow(slope_Y / (fbl_y * 8), 2));
                    slope = Math.Atan(slope) * 180 / Math.PI;
                    slope_res[(i - 1) + (j - 1) * dt.RasterXSize] = (int)(slope);
                    //aspect
                    double aspect = 180 / Math.PI * Math.Atan2(slope_Y/8, -slope_X/8);
                    if (slope_X == 0 && slope_Y == 0)//如果两个方向变化率为0，坡向赋值为-1
                    {
                        aspect = -1;
                    }
                    else if (aspect < 0)
                    {
                        aspect = 90.0 - aspect;

                    }else if (aspect > 90)
                    {
                        aspect = 450.0 - aspect;

                    }else 
                    { 
                        aspect = 90.0 - aspect;
                    }
                    aspect_res[(i - 1) + (j - 1) * dt.RasterXSize] = (int)(aspect);
                }
            }
        }
        //计算坡度坡向
        public void slopeAndAspect(Dataset dt,string path1,string path2)
        {
            // 获取仿射变换六参数
            double[] trans = new double[6];
            dt.GetGeoTransform(trans);

            //获取横向、纵向分辨率
            double fbl_x = trans[1];
            double fbl_y = trans[5];

            //定义模板(分别定义x方向和y方向的坡度模板)
            int[,] ope_x = new int[3, 3];
            ope_x[0, 0] = -1; ope_x[0, 1] = 0; ope_x[0, 2] = 1;
            ope_x[1, 0] = -2; ope_x[1, 1] = 0; ope_x[1, 2] = 2;
            ope_x[2, 0] = -1; ope_x[2, 1] = 0; ope_x[2, 2] = 1;

            int[,] ope_y= new int[3,3];
            ope_y[0,0] = -1; ope_y[0,1] =-2;  ope_y[0,2] =-1;
            ope_y[1,0] =  0; ope_y[1,1] = 0;  ope_y[1,2] =0;
            ope_y[2,0] =  1; ope_y[2,1] = 2;  ope_y[2,2] =1;

            //读取dem数据并转化为二维数组
            double[] dem = fulldem(dt);
            double[,] dem2D = To2D(dem, dt.RasterXSize+2, dt.RasterYSize+2);
            
            //求坡度坡向
            OperatorEnforce(dt,dem2D,ope_x,ope_y,fbl_x,fbl_y,out int[]slope,out int[]aspect);
            cls_saveFiles sf = new cls_saveFiles();

            List<int[]> slopes = new List<int[]>();
            List<int[]> aspects = new List<int[]>();

            slopes.Add(slope);
            aspects.Add(aspect);

            sf.SaveFromDataset(dt, path1, slopes, false);
            sf.SaveFromDataset(dt, path2, aspects, false);
        }
        public void slopeAndAspect2(Dataset dt, out int[] slope, out int[] aspect)
        {
            // 获取仿射变换六参数
            double[] trans = new double[6];
            dt.GetGeoTransform(trans);

            //获取横向、纵向分辨率
            double fbl_x = trans[1];
            double fbl_y = trans[5];

            //定义模板(分别定义x方向和y方向的坡度模板)
            int[,] ope_x = new int[3, 3];
            ope_x[0, 0] = -1; ope_x[0, 1] = 0; ope_x[0, 2] = 1;
            ope_x[1, 0] = -2; ope_x[1, 1] = 0; ope_x[1, 2] = 2;
            ope_x[2, 0] = -1; ope_x[2, 1] = 0; ope_x[2, 2] = 1;

            int[,] ope_y = new int[3, 3];
            ope_y[0, 0] = -1; ope_y[0, 1] = -2; ope_y[0, 2] = -1;
            ope_y[1, 0] = 0; ope_y[1, 1] = 0; ope_y[1, 2] = 0;
            ope_y[2, 0] = 1; ope_y[2, 1] = 2; ope_y[2, 2] = 1;

            //读取dem数据并转化为二维数组
            double[] dem = fulldem(dt);
            double[,] dem2D = To2D(dem, dt.RasterXSize + 2, dt.RasterYSize + 2);

            //求坡度坡向
            OperatorEnforce(dt, dem2D, ope_x, ope_y, fbl_x, fbl_y,  out slope, out aspect);
        }
        //一维数组转为二维数组
        public double[,] To2D(double[] array1D, int xsize, int ysize)
        {
            double[,] array2D = new double[xsize, ysize];
            for (int i = 0; i < xsize; i++)//列
            {
                for (int j = 0; j < ysize; j++)//行
                {
                    array2D[i, j] = array1D[i + j * xsize];
                }
            }
            return array2D;
        }
        #endregion
       
        #region 余弦校正
        public List<int[]> cosCorrection(Dataset dem_dt, Dataset dt, double gdj, double fwj)
        {
            //事先要dem重采样（dem数据和原图的分辨率可能不同）
            //获取坡度坡向
            slopeAndAspect2(dem_dt, out int[] slope, out int[] aspect);

            //光线入射角
            double[] cos_rushejiao = new double[dem_dt.RasterXSize * dem_dt.RasterYSize];
            double zh = Math.PI / 180;

            //太阳高度角的余弦值
            double cos_gdj = Math.Cos(gdj * zh);
            double sin_gdj = Math.Sin(gdj * zh);

            //计算光线入射角的余弦值
            for (int i = 0; i < dem_dt.RasterXSize * dem_dt.RasterYSize; i++)
            {
                cos_rushejiao[i] = cos_gdj * Math.Cos(slope[i] * zh) + sin_gdj * Math.Cos((fwj - aspect[i]) * zh);
            }

            List<int[]> newvalue = new List<int[]>();
            cls_basicfunc bc = new cls_basicfunc();
            int bandnum = dt.RasterCount;
            int[,] old = new int[bandnum,dt.RasterXSize*dt.RasterYSize];
            List<int[]> oldvalue = bc.getvalue(dt);
            for(int count = 0; count < oldvalue.Count; count++)
            {
                for (int i = 0; i < dt.RasterXSize * dt.RasterYSize; i++)
                {
                    old[count, i] = oldvalue[count][i];
                }
            }
         
            //余弦校正
            for (int count = 0; count < oldvalue.Count; count++)
            {
                int[] newbuf = new int[dt.RasterXSize * dt.RasterYSize];
                for (int i = 0; i < dt.RasterXSize * dt.RasterYSize; i++)
                {
                    newbuf[i] = (int)(old[count,i] * cos_gdj / cos_rushejiao[i]);
                }
                newvalue.Add(newbuf);
            }
            return newvalue;

        }
        #endregion

        #region C校正
        public List<int[]> Ccorrection(Dataset dem_dt, Dataset dt, double gdj, double fwj,out double []a,out double []b)
        {
            //事先要dem重采样（dem数据和原图的分辨率可能不同）
            //获取坡度坡向
            slopeAndAspect2(dem_dt, out int[] slope, out int[] aspect);

            //光线入射角
            double[] cos_rushejiao = new double[dem_dt.RasterXSize * dem_dt.RasterYSize];
            double zh = Math.PI / 180;

            //太阳高度角的余弦值
            double cos_gdj = Math.Cos(gdj * zh);
            double sin_gdj = Math.Sin(gdj * zh);

            //计算光线入射角的余弦值
            for (int i = 0; i < dem_dt.RasterXSize * dem_dt.RasterYSize; i++)
            {
                cos_rushejiao[i] = cos_gdj * Math.Cos(slope[i] * zh) + sin_gdj * Math.Cos((fwj - aspect[i]) * zh);
            }

            List<int[]> newvalue = new List<int[]>();
            cls_basicfunc bc = new cls_basicfunc();

            int[,] old = new int[dt.RasterCount, dt.RasterXSize * dt.RasterYSize];
            List<int[]> oldvalue = bc.getvalue(dt);

            for (int count = 0; count < oldvalue.Count; count++)
            {
                for (int i = 0; i < dt.RasterXSize * dt.RasterYSize; i++)
                {
                    old[count, i] = oldvalue[count][i];
                }
            }

            //确定拟合的系数(最小二乘线性拟合)

            cls_basicfunc bf = new cls_basicfunc();
            var cos_i = new DenseMatrix(cos_rushejiao.Length, 1);
            var Lt = new DenseMatrix(cos_rushejiao.Length, 1);
            List<Matrix> css = new List<Matrix>();
            //取样本：取中间位置的一行、一列
            for (int i = 0; i < dt.RasterXSize; i++)
            {
                cos_i[i,0]= cos_rushejiao[i+dt.RasterXSize*(dt.RasterYSize/2)];
            }
            for (int j = 0; j < dt.RasterYSize; j++)
            {
                cos_i[j+dt.RasterXSize, 0] = cos_rushejiao[j*dt.RasterXSize + dt.RasterXSize / 2];
            }

            //取样本：取中间位置的一行、一列
            for (int num = 0; num < dt.RasterCount;num++)
            {
                for (int i = 0; i < dt.RasterXSize; i++)
                {
                    Lt[i, 0] = oldvalue[num][i + dt.RasterXSize * (dt.RasterYSize / 2)];
                }
                for (int j = 0; j < dt.RasterYSize; j++)
                {
                    Lt[j + dt.RasterXSize, 0] = oldvalue[num][j * dt.RasterXSize + dt.RasterXSize / 2];
                }
                var cs = bf.linearRegression(cos_i,Lt);
                css.Add((Matrix)cs);
            }

            a = new double[dt.RasterCount];
            b = new double[dt.RasterCount];

            for (int num = 0; num < dt.RasterCount; num++)
            {
                a[num] = css[num].At(1,0);
                b[num] = css[num].At(0,0);
            }

           //C校正
                for (int count = 0; count < oldvalue.Count; count++)
            {
                int[] newbuf = new int[dt.RasterXSize * dt.RasterYSize];
                for (int i = 0; i < dt.RasterXSize * dt.RasterYSize; i++)
                {
                    newbuf[i] = (int)(old[count, i] * (cos_gdj+a[count]/b[count]) /(cos_rushejiao[i]+a[count]/b[count]));
                }
                newvalue.Add(newbuf);
            }
            return newvalue;

        }
        #endregion
     
    }

}
