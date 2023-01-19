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

        #region GMM 高斯混合模型(1) 待更正
        //public double[] initialization_pi(int class_num)//初始化权重pi矩阵
        //{
        //    double[] pi = new double[class_num];
        //    Random ran = new Random();
        //    double sum = 0.0;
        //    for (int i=0;i<class_num;i++)//随机数
        //    {
        //        pi[i] = ran.Next(0,100);
        //        sum += pi[i];
        //    }
        //    for (int i = 0; i < class_num; i++)//概率归一化
        //    {
        //        pi[i] /=sum;
        //    }
        //    return pi;
        //}
        //public double[,] initialization_miu(Dataset dt, int class_num)//初始化均值miu矩阵
        //{
        //    int bands = dt.RasterCount;
        //    double[,] miu = new double[class_num,bands];
        //    Random ran = new Random();

        //    //求图像各波段像素值范围
        //    cls_basicfunc bc = new cls_basicfunc();
        //    int[] min = bc.MaxMin(dt, "min");
        //    int[] max = bc.MaxMin(dt, "max");
        //    int[] range = new int[bands];
        //    for(int i=0;i<bands;i++)
        //    {
        //        range[i] = max[i] - min[i];
        //    }

        //    //分配不同类的均值
        //    for (int Class=0;Class<class_num;Class++)
        //    {
        //        for(int band=0;band<bands;band++)
        //        {
        //            miu[Class, band] = min[band] + range[band] * ran.NextDouble();
        //        }
        //    }
        //    return miu;
        //}
        //public Matrix<double> toMiuk(double[,] miu, int classflag, int bandnum)// 求出miu（k），即第k类的均值矩阵
        //{
        //    var miuk = new DenseMatrix(bandnum, 1);
        //    for (int i = 0; i < bandnum; i++)
        //    {
        //        miuk[i, 0] = miu[classflag, i];
        //    }
        //    return miuk;
        //}
        //public double[,,]initialization_sigma(Dataset dt, int class_num)//初始化方差sigma矩阵
        //{
        //    int bands = dt.RasterCount;
        //    double[,,] sigma = new double[class_num,bands,bands];
        //    Random ran = new Random();

        //    //求图像各波段间协方差
        //    cls_basicfunc bc = new cls_basicfunc();
        //    double [,]xfc =bc.xfc(dt);

        //    //初始化不同类的协方差
        //    for (int Class = 0; Class < class_num; Class++)
        //    {
        //        for (int i = 0; i < bands; i++)
        //        {
        //            for(int j= 0; j < bands; j++)
        //            {
        //                sigma[Class,i,j] = xfc[i,j]*ran.NextDouble();
        //            }
        //        }
        //    }
        //    return sigma;
        //}
        //public Matrix<double> toSk(double[,,] sigma, int classflag, int bandnum)// 求出sigma(k),即第k类的协方差阵
        //{
        //    var sk = new DenseMatrix(bandnum, bandnum);
        //    for (int i = 0; i < bandnum; i++)
        //    {
        //        for (int j = 0; j < bandnum; j++)
        //        {
        //            sk[i, j] = sigma[classflag, i,j];
        //        }
        //    }
        //    return sk;
        //}
        //public Matrix<double> toXi(int [,]X , int pixelflag,int bandnum)// 求出Xi(i像素的所有band组成的向量)
        //{
        //    var x = new DenseMatrix(bandnum,1);
        //    for(int i=0;i<bandnum;i++)
        //    {
        //        x[i, 0] = X[pixelflag, i];
        //    }
        //    return x;
        //}
        //public double Gauss(Matrix<double>Xi, Matrix<double> miuk,Matrix<double>sigmak) //多维高斯概率密度分布函数
        //{
        //    double N=0.0;
        //    double a = 0.0, aa=0.0, aaa=0.0;
        //    double num_zs = 0.0;

        //    double d = (double)(Xi.RowCount);

        //    double b = sigmak.Determinant();
        //    aa = Math.Sqrt(Math.Abs(b));
        //    aaa = Math.Pow((2 * Math.PI), d / 2);
        //    a = 1 / aaa / aa;

        //    var zs = new DenseMatrix(1,1);
        //    zs =(DenseMatrix)( -0.5 * (Xi - miuk).Transpose() * sigmak.Inverse() * (Xi - miuk));
        //    num_zs = zs[0, 0];

        //    N = a * Math.Exp(num_zs);

        //    return N;
        //}
        //public double[,] posterP(Dataset dt, int[,]X, double[,]miu,double[,,]sigma,double [] pi,int class_num)//计算后验概率矩阵P
        //{
        //    int pixels = dt.RasterXSize * dt.RasterYSize;
        //    int bandnum = dt.RasterCount;
        //    double[,] poster = new double[pixels,class_num];

        //    for(int p=0;p<pixels;p++)
        //    {
        //        double N = 0.0;
        //        double sum = 0.0;
        //        Matrix<double> Xi = toXi(X,p,bandnum);

        //        for (int c=0;c<class_num;c++)
        //        {
        //            Matrix<double> miuk = toMiuk(miu, c, bandnum);
        //            Matrix<double> simgak = toSk(sigma, c, bandnum);
        //            N = pi[c] * Gauss(Xi,miuk,simgak);

        //            for (int cc = 0; cc < class_num; cc++)
        //            {
        //                Matrix<double> miukk = toMiuk(miu, c, bandnum);
        //                Matrix<double> simgakk = toSk(sigma, c, bandnum);
        //                N = pi[c] * Gauss(Xi, miukk, simgakk);
        //                sum += N;
        //            }

        //            poster[p, c] = N / sum;
        //        }
        //    }

        //    return poster;

        //}
        //public double[,] new_miu(Dataset dt,int[,]X,double[,]p,int class_num)//更新均值矩阵
        //{
        //    int pixels = dt.RasterXSize * dt.RasterYSize;
        //    int bands = dt.RasterCount;
        //    double[,] new_miu = new double[class_num,bands];

        //    for (int k = 0; k < class_num; k++)
        //    {
        //        double Nk = 0.0;
        //        for (int i = 0; i < pixels; i++)
        //        {
        //            Nk += p[i, k];
        //        }

        //        for (int b = 0; b < bands; b++)
        //        {
        //            double sum = 0.0;
        //            for (int i = 0; i < pixels; i++)
        //            {
        //                 sum +=  p[i, k] * X[i,b];
        //            }
        //            new_miu[k, b] = 1 / Nk * sum;
        //        }
        //    }
        //    return new_miu;
        //}
        //public double[] new_pi(Dataset dt, double[,] p, int class_num)//更新权重矩阵
        //{
        //    int pixels = dt.RasterXSize * dt.RasterYSize;
        //    double[] new_pi = new double[class_num];

        //    for (int k = 0; k < class_num; k++)
        //    {
        //        double Nk = 0.0;
        //        for (int i = 0; i < pixels; i++)
        //        {
        //            Nk += p[i, k];
        //        }
        //        new_pi[k] = Nk / pixels;
        //    }
        //    return new_pi;
        //}
        //public double[,,] new_sigma(Dataset dt, int[,] X, double[,] p, double[,]miu, int class_num)//更新协方差矩阵
        //{
        //    int pixels = dt.RasterXSize * dt.RasterYSize;
        //    int bands = dt.RasterCount;
        //    double[,,] new_sigma = new double[class_num,bands,bands];

        //    for(int k=0;k<class_num;k++)
        //    {
        //        double Nk = 0.0;
        //        var miuk = toMiuk(miu, k, bands);
        //        var sigmak = DenseMatrix.CreateDiagonal(bands,bands,0.0);

        //        for (int i = 0; i < pixels; i++)//求Nk
        //        {
        //            Nk += p[i, k];
        //        }

        //        for(int j=0;j<pixels;j++)//求sigma(k)
        //        {
        //            var xi = toXi(X,j,bands);
        //            sigmak += (DenseMatrix)(1 / Nk * p[j, k] * (xi-miuk) * (xi - miuk).Transpose());
        //        }

        //        for(int i=0;i<bands;i++)//求sigma
        //        {
        //            for(int j=0;j<bands;j++)
        //            {
        //                new_sigma[k, i, j] = sigmak[i,j];
        //            }
        //        }
        //    }
        //    return new_sigma;
        //}
        //public void GMM(Dataset dt,int class_num,int loop,string filepath)  // 混合高斯模型（GMM） 
        //{
        //    cls_basicfunc bc = new cls_basicfunc();
        //    List<int[]> data = bc.getvalue(dt);

        //    int pixels = data[0].Length;
        //    int bands = data.Count;

        //    // 构建光谱特征向量 [像元，波段号]
        //    int[,] X = new int[pixels, bands];
        //    for (int i = 0; i < pixels; i++)
        //    {
        //        for (int j = 0; j < bands; j++)
        //        {
        //            X[i, j] = data[j][i];
        //        }
        //    }

        //    //初始化参量
        //    double[] pi = initialization_pi(class_num);
        //    double[,] miu = initialization_miu(dt, class_num);
        //    double[,,] sigma = initialization_sigma(dt,class_num);

        //    //迭代
        //    double[,] p = new double[pixels, class_num];
        //    for(int i=0;i<loop;i++)
        //    {
        //        //计算后验概率p
        //        p = posterP(dt,X,miu,sigma,pi,class_num);
        //        // 计算新的参量
        //        miu = new_miu(dt, X, p, class_num);
        //        sigma = new_sigma(dt, X, p, miu, class_num);
        //        pi = new_pi(dt, p, class_num);

        //        double sum = 0.0;
        //        for(int j=0;j<class_num;j++)
        //        {
        //            sum += pi[j];
        //        }
        //        for (int j = 0; j < class_num; j++)
        //        {
        //            pi[j]/= sum ;
        //        }
        //    }

        //    //获取每个像元所属的类别
        //    int[] classes = new int[pixels];
        //    for (int j = 0; j < pixels; j++)
        //    {
        //        int classflag = 0;
        //        double max = p[j, 0];
        //        for (int c = 1; c < class_num; c++)
        //        {
        //            if (p[j, c] > max)
        //            {
        //                max = p[j, c];
        //                classflag = c;
        //            }
        //        }
        //        classes[j] = classflag;
        //    }

        //    //给分类后图像赋色(灰度图)
        //    int[] value = new int[pixels];//赋值后的每个像素值的集合
        //    if (class_num == 2)
        //    {
        //        for (int j = 0; j < pixels; j++)
        //        {
        //            if (classes[j] == 0)
        //                value[j] = 0;
        //            if (classes[j] == 1)
        //                value[j] = 255;
        //        }
        //    }

        //    int[] temp1 = new int[class_num];
        //    temp1[0] = 0;
        //    temp1[class_num - 1] = 255;
        //    if (class_num > 2)
        //    {
        //        for (int i = 1; i < class_num - 1; i++)
        //        {
        //            temp1[i] = temp1[i - 1] + 256 / (class_num - 1);
        //        }

        //        for (int j = 0; j < pixels; j++)
        //        {
        //            value[j] = temp1[classes[j]];
        //        }
        //    }

        //    cls_saveFiles sf = new cls_saveFiles();
        //    List<int[]> dns = new List<int[]> { value };
        //    sf.SaveFromDataset(dt, filepath, dns, false);
        //}
        #endregion

        #region GMM 高斯混合模型(2) 可用

        //初始化均值矩阵
        public double[,] Initialization_miu(Dataset dt, int K)
        {
            int bands = dt.RasterCount;
            double[,] u = new double[K, bands];
            Random ran = new Random();

            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < bands; j++)
                {
                    Band band = dt.GetRasterBand(j + 1);
                    double[] maxandmin = { 0, 0 };
                    band.ComputeRasterMinMax(maxandmin, 0);
                    double randNum = ran.NextDouble();
                    u[i, j] = maxandmin[0] + (maxandmin[1] - maxandmin[0]) * randNum;
                }
            }
            return u;
        }
        
        //计算像元xi属于第k类的高斯概率
        public double pro(double[] x, int k, double[,] m_vars, double[,] u)
        {
            double p = 1;
            int dimNum = x.Length;
            for (int d = 0; d < dimNum; d++)
            {
                p *= 1 / Math.Sqrt(2 * Math.PI * m_vars[k, d]);
                p *= Math.Exp(-0.5 * (x[d] - u[k, d]) * (x[d] - u[k, d]) / m_vars[k, d]);
            }
            return p;
        }

        //计算像元xi属于所有类的总高斯概率
        public double totalpro(double[] x, int K, double[] weights, double[,] m_vars, double[,] u)
        {
            double p = 0;
            for (int i = 0; i < K; i++)
            {
                p += weights[i] * pro(x, i, m_vars, u);
            }
            return p;
        }
        
        public void GMM2(Dataset dt, int K, int loop, string filepath)
        {
            re:
            cls_basicfunc bc = new cls_basicfunc();
            List<int[]> data = bc.getvalue(dt);//获取原数据DN值
            int pixels = data[0].Length;//像素个数
            int bands = dt.RasterCount;//波段数

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
            double[,] u = Initialization_miu(dt, K);//均值初始化
            int[] types = new int[pixels];

            //根据均值对像元进行初始分类
            for (int i = 0; i < pixels; i++)
            {
                double mindistance = dist(X, u, bands, i, 0);
                for (int j = 0; j < K; j++)
                {
                    if (mindistance > dist(X, u, bands, i, j))
                    {
                        mindistance = dist(X, u, bands, i, j);
                        types[i] = j;
                    }
                }
            }

            //每一类的数量
            double[] counts = new double[K];
            for (int i = 0; i < types.Length; i++)
            {
                counts[types[i]]++;
            }

            // 计算先验概率权重
            double[] weights = new double[K];
            for (int i = 0; i < K; i++)
            {
                weights[i] = counts[i] / pixels;
            }

            double[,] vars = new double[K, bands];
            // 计算每个分类的方差
            double[] overMeans = new double[bands];
            double[] m_minVars = new double[bands];
            double[] x = new double[bands];
            for (int i = 0; i < pixels; i++)
            {
                for (int j = 0; j < bands; j++)
                {
                    x[j] = X[i, j];
                }

                for (int d = 0; d < bands; d++)
                {
                    vars[types[i], d] += (x[d] - u[types[i], d]) * (x[d] - u[types[i], d]);
                }
             
                //计算总体均值和方差
                for (int d = 0; d < bands; d++)
                {
                    overMeans[d] += x[d];
                    m_minVars[d] += x[d] * x[d];
                }
            }

            double MIN_VAR = 1E-10;
            // Compute the overall variance (* 0.01) as the minimum variance.
            for (int d = 0; d < bands; d++)
            {
                overMeans[d] /= pixels;
                m_minVars[d] = Math.Max(MIN_VAR, 0.01 * (m_minVars[d] / pixels - overMeans[d] * overMeans[d]));
            }

            // Initialize each Gaussian.
            for (int i = 0; i < K; i++)
            {
                if (weights[i] > 0)
                {
                    for (int d = 0; d < bands; d++)
                    {
                        vars[i, d] = vars[i, d] / counts[i];

                        // A minimum variance for each dimension is required.
                        if (vars[i, d] < m_minVars[d])
                        {
                            vars[i, d] = m_minVars[d];
                        }
                    }
                }
            }

            double[,] next_means = new double[K, bands];
            double[] next_weights = new double[K];
            double[,] next_vars = new double[K, bands];
            int[] pixelofclass = new int[pixels];//记录像元类别
            int times = 0;

            //进行迭代运算
            while (times < loop)
            {
                //首先初始化要用的数组
                for (int j1 = 0; j1 < K; j1++)
                {
                    next_weights[j1] = 0;
                    for (int j2 = 0; j2 < bands; j2++)
                    {
                        next_means[j1, j2] = 0;
                        next_vars[j1, j2] = 0;
                    }
                }

                for (int k = 0; k < pixels; k++)
                {
                    for (int j = 0; j < bands; j++)
                        x[j] = X[k, j];

                    double p = totalpro(x, K, weights, vars, u); // 总的概率密度分布
                    double maxp = 0;
                    for (int j = 0; j < K; j++)
                    {
                        double pj = pro(x, j, vars, u) * weights[j] / p; // 每个分类的概率密度分布百分比
                        if (maxp < pj)
                        {
                            maxp = pj;
                            pixelofclass[k] = j;
                        }

                        next_weights[j] += pj; // 得到后验概率

                        for (int d = 0; d < bands; d++)
                        {
                            next_means[j, d] += pj * x[d];
                            next_vars[j, d] += pj * x[d] * x[d];
                        }
                    }

                }

                // 重新估计：生成新的权重、均值和方差
                for (int j = 0; j < K; j++)
                {
                    weights[j] = next_weights[j] / pixels;

                    if (weights[j] > 0)
                    {
                        for (int d = 0; d < bands; d++)
                        {
                            u[j, d] = next_means[j, d] / next_weights[j];
                            vars[j, d] = next_vars[j, d] / next_weights[j] - u[j, d] * u[j, d];
                            if (vars[j, d] < m_minVars[d])
                            {
                                vars[j, d] = m_minVars[d];
                            }
                        }
                    }
                }
                times++;//次数加一
            }

            // 输出限制：如果有类别的像元数量为0，重新分类
            int []count = new int [K];
            for(int i=0;i<pixels;i++)
            {
                count[pixelofclass[i]]++;
            }
            for(int i=0;i<K;i++)
            {
                if(count[i]==0)
                {
                    goto re;
                }
            }

            #region 输出灰度分类图
            //输出图像
            //给分类后图像赋色
            //灰度图
            int[] colors = new int[pixels];//赋值后的每个像素值的集合
            if (K == 2)
            {
                for (int j = 0; j < pixels; j++)
                {
                    if (pixelofclass[j] == 0)
                        colors[j] = 0;
                    if (pixelofclass[j] == 1)
                        colors[j] = 255;
                }
            }

            int[] temp1 = new int[K];
            temp1[0] = 0;
            temp1[K - 1] = 255;
            if (K > 2)//类别数大于2
            {
                for (int i = 1; i < K - 1; i++)
                {
                    temp1[i] = temp1[i - 1] + 256 / (K - 1);
                }

                for (int j = 0; j < pixels; j++)
                {
                    colors[j] = temp1[pixelofclass[j]];
                }
            }

            List<int[]> dns = new List<int[]> { colors };
            cls_saveFiles sf = new cls_saveFiles();
            sf.SaveFromDataset(dt, filepath, dns, false);
            #endregion
        }

        #endregion

        #region 最小距离分类
        public void mindisClassify(Dataset dt, double[,]miu,int classnum , List<Color>colors,string filepath)
        {
            cls_basicfunc bc = new cls_basicfunc();
            List<int[]> data = bc.getvalue(dt);
            int bandnum = data.Count;//波段数
            int pixels = data[0].Length;//像素数
            var Xi = new DenseMatrix(bandnum,1);//第i个像元的特征向量
            var miuk = new DenseMatrix(bandnum,1);//第k类的均值向量
            double[] dis = new double[classnum];//距离
            int[] pixels_of_class = new int[pixels];//每个像元所属类别
            //找出每一个像元所属的类别
            for(int i=0;i<pixels;i++)
            {
                for(int j=0;j<bandnum;j++)
                {
                    Xi[j, 0] = data[i][j];
                }

                for(int k=0;k<classnum;k++)
                {
                    for(int j=0;j<bandnum;j++)
                    {
                        miuk[j, 0] = miu[k,j];
                    }
                    var diss = (Xi - miuk).Transpose() * (Xi - miuk);
                    dis[k] = diss[0,0];
                }
                double min = dis[0];
                int c = 0;
                for(int k=0;k<classnum;k++)
                {
                    if (dis[k] < min)
                    {
                        min = dis[k];
                        c = k;
                    }
                }
                pixels_of_class[i] = c;
            }
            List<int[]> RGB = new List<int[]>();
            for(int i=0;i<3;i++)
            {
                RGB.Add(new int [pixels]);
            }
            //赋色
            for(int i=0;i<pixels;i++)
            {
                RGB[0][i] = colors[pixels_of_class[i]].R;
                RGB[1][i] = colors[pixels_of_class[i]].G;
                RGB[2][i] = colors[pixels_of_class[i]].B;
            }
            //保存
            cls_saveFiles sf = new cls_saveFiles();
            sf.SaveFromDataset(dt, filepath, RGB, false);

        }
        #endregion

    }
}
