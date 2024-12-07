# seven-sorting-ways
#十大排序之七（不含冒泡，插入，选择）  
我在cpp文件里面写了7种排序方法，并通过代码计算了它们的性能差异  

具体排名如下：
综合排名(10数量级权重0.05,100数量级权重0.15，高数量级权重0.8）  
No.1 基数( 3*0.05 + 7*0.15 + 6*0.8 = 6.00 )  
No.2 计数( 1*0.05 + 2*0.15 + 7*0.8 = 5.95 )  
No.3 快速( 7*0.05 + 6*0.15 + 5*0.8 = 5.25 )  
No.4 归并( 5*0.05 + 5*0.15 + 4*0.8 = 4.20 )  
No.5 堆  ( 6*0.05 + 4*0.15 + 3*0.8 = 3.30 )  
No.6 希尔( 4*0.05 + 3*0.15 + 2*0.8 = 2.25 )  
No.7 桶  ( 2*0.05 + 1*0.15 + 1*0.8 = 1.05 )  

你可以查看我的排序代码，通过修改排序的数组长度和排序范围来查看不同排序方式的时间  
