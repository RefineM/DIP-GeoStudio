using OSGeo.GDAL;
using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;
using System.Drawing;

namespace dip_1
{
    class cls_classify
    {

        #region FCM 模糊C均值分类
        public double dist(int[,] X, double[,] C, int bands, int j, int i)//求向量欧氏距离
        {
            double dis = 0.0;
            double temp = 0.0;
            for (int band = 0; band < bands; band++)
            {
                temp += Math.Pow((X[j, band] - C[i, band]), 2);
            }
            dis = Math.Sqrt(temp);
            return dis;
        }
        public double[,] initializtion_U(int pixels, int class_num)//初始化隶属度矩阵
        {
            double[,] U = new double[pixels, class_num];
            double[] sum = new double[pixels];

            //生成随机数
            Random ran = new Random();
            for (int i = 0; i < pixels; i++)
                for (int j = 0; j < class_num; j++)
                {
                    U[i, j] = ran.Next(0, 100);
                }

            //求出每一行的和
            for (int i = 0; i < pixels; i++)
                for (int j = 0; j < class_num; j++)
                {
                    sum[i] += U[i,j];
                }

            //把概率归一化到【0，1】区间
            for (int i = 0; i < pixels; i++)
                for (int j = 0; j < class_num; j++)
                {
                    U[i, j] = U[i, j] / sum[i];
                }
            return U;
        }
        public void FCM(Dataset dt, double m, int class_num, int loop,string filepath)// 模糊C均值聚类（FCM）
        {
            cls_basicfunc bc = new cls_basicfunc();
            List<int[]> data = bc.getvalue(dt);
            int pixels = data[0].Length;
            int bands = data.Count;

            // 构建光谱特征向量 [像元，波段号]
            int[,] X = new int[pixels, bands];
            //获取向量
            for (int i = 0; i < pixels; i++)
            {
                for (int j = 0; j < bands; j++)
                {
                    X[i, j] = data[j][i];
                }
            }
            //初始化隶属度矩阵
            double[,] U = initializtion_U(pixels, class_num);
            //聚类中心矩阵
            double[,] C = new double[class_num, bands];

            //迭代
            for (int k = 0; k <= loop; k++)
            { 

                //计算聚类中心 C（i，b） [class,bands]
                for (int band = 0; band < bands; band++)
                {
                    for (int i = 0; i < class_num; i++)
                    {
                        double u_sum = 0.0;
                        double xu_sum = 0.0;

                        for (int j = 0; j < pixels; j++)
                        {
                            u_sum += Math.Pow(U[j, i], m);
                        }

                        for (int j = 0; j < pixels; j++)
                        {
                            xu_sum += Math.Pow(U[j, i], m) * X[j, band];
                        }

                        C[i, band] = xu_sum / u_sum;
                    }    
                }

                ////计算目标函数 J 
                //double dis = 0.0;
                //for (int i = 0; i < class_num; i++)
                //{
                //    for (int j = 0; j < pixels; j++)
                //    {
                //        double temp = 0.0;
                //        for (int band = 0; band < bands; band++)
                //        {
                //            temp += Math.Pow((X[j, band] - C[i, band]), 2);
                //        }
                //        dis = temp * U[j, i];
                //        J[k] += dis;
                //    }
                //}

                //更新隶属度矩阵U
                for (int i = 0; i < class_num; i++)
                {
                    for (int j = 0; j < pixels; j++)
                    {
                        double trans = 0.0;
                        for (int kk = 0; kk < class_num; kk++)
                        {
                            trans += Math.Pow((dist(X, C, bands, j, i) / dist(X, C, bands, j, kk)), (2 / (m - 1)));
                        }
                        U[j, i] = 1 / trans;
                    }
                }
            }

            //获取每个像元所属的类别
            int[] classes = new int[pixels];
            for (int j = 0; j < pixels; j++)
            {
                int classflag = 0;
                double max = U[j, 0];
                for (int c = 1; c < class_num; c++)
                {
                    if (U[j, c] > max)
                    {
                        max = U[j, c];
                        classflag = c;
                    }
                }
                classes[j] = classflag;
            }

            //给分类后图像赋色
            //灰度图
            int[] value = new int[pixels];//赋值后的每个像素值的集合
            if (class_num == 2)
            {
                for (int j = 0; j < pixels; j++)
                {
                    if (classes[j] == 0)
                        value[j] = 0;
                    if (classes[j] == 1)
                        value[j] = 255;
                }
            }
            
            int[] temp1 = new int[class_num];
            temp1[0] = 0;
            temp1[class_num - 1] = 255;
            if (class_num > 2)
            {
                for (int i = 1; i < class_num - 1; i++)
                {
                    temp1[i] = temp1[i - 1] + 256 / (class_num - 1);
                }

                for (int j = 0; j < pixels; j++)
                {
                    value[j] = temp1[classes[j]];
                }
            }

            cls_saveFiles sf = new cls_saveFiles();
            List<int[]> dns = new List<int[]> {value};
            sf.SaveFromDataset(dt, filepath, dns, false);
        }
        #endregion

        #region K-means
        public void kmeans (Dataset dt,int classnum,string path)
        {
            cls_basicfunc bc = new cls_basicfunc();
            List<int[]> data = bc.getvalue(dt);
            int bands = data.Count;//波段数
            int pixels = data[0].Length;//像素数

            // 构建光谱特征向量 [像元，波段号]
            int[,] X = new int[pixels, bands];
            // 获取向量
            for (int i = 0; i < pixels; i++)
            {
                for (int j = 0; j < bands; j++)
                {
                    X[i, j] = data[j][i];
                }
            }
            re:
            // 建立聚类中心并初始化
            double[,]C = initial_c(dt,bands,classnum);
            //每个像素所属类别号存储
            int[] pixel_class = new int[pixels];
          
            while(1==1)
            {
                //计算当前每个像元所属类别
                for (int i=0;i<pixels;i++)
                {
                    double[] dis = new double[classnum];
                    for (int k=0;k<classnum;k++)
                    {
                        dis[k] = distance(X,C,bands,i,k);
                    }

                    double min = dis[0];
                    int flag = 0;
                    for(int k = 0; k < classnum; k++)
                    {
                        if(dis[k]<min)
                        {
                            min = dis[k];
                            flag = k;
                        }
                    }

                    pixel_class[i] = flag;
                }

                //每一类的像元数,总和
                int[] pixel_num = new int[classnum];
                int[,] sum = new int[classnum,bands];

                for (int i = 0; i < pixels; i++)
                {
                    pixel_num[pixel_class[i]]++;
                }

                for(int b = 0; b < classnum; b++)
                {
                    if(pixel_num[b]==0)
                    {
                        goto re;
                    }
                }

                for (int i = 0; i < pixels; i++)
                {
                    for (int b = 0; b < bands; b++)
                    {
                        sum[pixel_class[i], b] += X[i,b];
                    }
                }

                double[,] C_new = new double[classnum,bands];
                //新的聚类中心
                for(int c=0;c<classnum;c++)
                {
                    for(int b=0;b<bands;b++)
                    {
                        C_new[c, b] = sum[c, b] / pixel_num[c];
                    }
                }

                //计算前后两个聚类中心的距离
                double[] Cdis = new double[classnum];
                for(int i=0;i<classnum;i++)
                {
                     Cdis[i] = Cdistance(C,C_new,bands,i,i);
                }

                if(cls_fusion.maxvalueD(Cdis)<= 0.01)//迭代终止条件
                {
                    break;
                }
                else
                {
                    C = C_new;
                }

            }

            //输出灰度图像
            //根据类别，生成颜色表
            int[]color = new int[pixels];
            if (classnum == 2)
            {
                for (int j = 0; j < pixels; j++)
                {
                    if (pixel_class[j] == 0)
                        color[j] = 0;
                    if (pixel_class[j] == 1)
                        color[j] = 255;
                }
            }

            int[] temp1 = new int[classnum];
            temp1[0] = 0;
            temp1[classnum - 1] = 255;
            if (classnum > 2)
            {
                for (int i = 1; i < classnum - 1; i++)
                {
                    temp1[i] = temp1[i - 1] + 256 / (classnum - 1);
                }

                for (int j = 0; j < pixels; j++)
                {
                    color[j] = temp1[pixel_class[j]];
                }
            }

            cls_saveFiles sf = new cls_saveFiles();
            List<int[]> dns = new List<int[]> { color };
            sf.SaveFromDataset(dt, path, dns, false);

        }

        public double[,] initial_c(Dataset dt, int bands, int class_num)//初始化聚类
        {
            double[,] C = new double[class_num, bands];
            Random ran = new Random();

            //求图像各波段像素值范围
            cls_basicfunc bc = new cls_basicfunc();
            int[] min = bc.MaxMin(dt, "min");
            int[] max = bc.MaxMin(dt, "max");
            int[] range = new int[bands];
            for (int i = 0; i < bands; i++)
            {
                range[i] = max[i] - min[i];
            }

            //随机化
            for (int c=0;c<class_num;c++)
            {
                for(int i=0;i<bands;i++)
                {
                    C[c, i] = min[i] + range[i] * ran.NextDouble();
                }
            }
            return C;
        }

        public double distance(int[,]X,double[,]C,int bands,int xi,int ck )//计算向量间欧氏距离
        {
            double dis = 0.0;
            double temp = 0.0;
            for(int j=0;j<bands;j++)
            {
                temp += Math.Pow((X[xi,j]-C[ck,j]),2);
            }
            dis = Math.Sqrt(temp);
            return dis;
        }
        public double Cdistance(double[,] X, double[,] C, int bands, int xi, int ck)//计算向量间欧氏距离
        {
            double dis = 0.0;
            double temp = 0.0;
            for (int j = 0; j < bands; j++)
            {
                temp += Math.Pow((X[xi, j] - C[ck, j]), 2);
            }
            dis = Math.Sqrt(temp);
            return dis;
        }

        #endregion
      }
}
