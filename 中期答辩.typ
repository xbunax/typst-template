#import "@preview/polylux:0.3.1": *
#import "@preview/xarrow:0.2.0": xarrow, xarrowSquiggly, xarrowTwoHead
#import "@preview/whalogen:0.1.0": ce
#import "@preview/fletcher:0.2.0" as fletcher: node, edge
#import themes.university: *

#show: university-theme.with(short-author: "", short-title: "", short-date: "")
// #set text(size: 18pt)
//
#title-slide(
  title: "自旋轨道矩诱导的反铁磁磁矩翻转",
  // date: "2023年12月",
  authors: ("答辩人：朱墨" ),
  institution-name: text(font: "Zapfino", "Tongji University"),
)

#set text(font: "STFangsong")
#show figure.caption: it => [
  #set align(center)
  #text(size: 13pt, [#it.supplement #it.counter.display() #it.body])
]
//
// #slide(
//   title: [目录],
//   new-section: [],
// )[
//   #set text(size: 25pt, font: "STFangsong", weight: 500)
//   #let cell = rect.with(
//     width: 90%,
//     height: 15%,
//     stroke: none,
//     radius: (rest: 10pt),
//     fill: gray,
//   )
//   #let large_cell = rect.with(
//     width: 90%,
//     height: 40%,
//     stroke: none,
//     radius: (rest: 10pt),
//     fill: gray,
//   )
//   #set list(body-indent: 10pt)
//   #v(1.5em)
//   #grid(rows: 3, row-gutter: 10pt, cell[#align(horizon, [- 研究计划])], large_cell[
//     #set align(horizon)
//     - 自旋轨道矩诱导的反铁磁磁矩翻转
//     #set text(size: 20pt)
//     #list(indent: 25pt, marker: [], [
//       + DMI模拟程序
//       + 退磁场计算并行优化
//       + 程序验证
//     ])
//   ], cell[#align(horizon, [- 研究成果及下一步研究方向])])
//
// ]

// #slide(
//   title: [研究计划],
//   new-section: [],
// )[
//   #set align(horizon + center)
//   #let cell = rect.with(fill: gray, width: 95%, height: 49%, radius: (rest: 10pt))
//  #set text(size: 20pt)
//
//   #grid(
//     columns: 2,
//     rows: 2,
//     row-gutter: 5pt,
//     column-gutter: 0pt,
//     cell[#set align(left)
//       + *第一阶段(2022 年 9 月-2023 年 6 月)*
//
//       - 通过阅读书籍和文献，学习自旋电子学理论知识，掌握 STT 效应、SOT 效应原理，调研反铁磁相关研究的进展，并熟练掌握微 C 语言，有一定的编程能力。],
//     cell[#set align(left + top)
//       2. 第二阶段(2023年6月-2024年3月)
//
//       - 开始熟悉 `mmag2.0` 的程序结构，并阅读源码，了解 `mmag2.0` 如何运行。],
//     cell[#set align(left + top)
//       3. 第三阶段(2024年 3 月-2024年12月)
//
//       - 利用 CUDA 完成 DMI 模拟模块的编写，并利用现有的结果对模块进行 check，最后利用 DMI 模块进行科学计算。],
//     cell[#set align(left + top)
//       4. 第四阶段(2025 年 1 月-2025 年 3 月)
//
//       - 工作总结、撰写毕业论文。],
//   )
// ]

#slide(
  title: "自旋轨道矩诱导的反铁磁磁矩翻转",
  new-section: "微磁学模拟",
)[

  #let cell = rect.with(stroke: none, width: 100%, height: 100%)
  #set text(size: 24pt)
  #grid(
    columns: (60%, 40%),
    column-gutter: 5pt,
    cell[- 系统自由能

      $ E &=E_("exc") + E_d + E_k +E_("Zeeman") + E_("DMI") $
      $ E_("exc")&=J arrow(m)_i dot arrow(m)_j\
      E_d &= -mu_0 attach(limits(M), t: arrow) dot attach(limits(H), t: arrow, br: "d")\
      E_K &= K sum_(i,j,k)[1-(m_(i,j,k) dot k)^2] delta V \
      E_("Zeeman")&=-mu_0 attach(limits(M), t: arrow) dot attach(limits(H), t: arrow, br: "Zeeman")\
      E_("DM") &= D sum arrow(r)_(i,j) dot (arrow(m)_i times arrow(m)_j)
      $

    ],
    cell[
      #set align(center + horizon)
      #set text(size: 15pt)
      #figure(image("pic/different_spin.png"), caption: [有限差分法示意图])
    ],
  )

]

#slide(
  title: "自旋轨道矩诱导的反铁磁磁矩翻转",
  new-section: "DMI模拟程序",
)[
  #set text(size: 18pt)
  #let cell = rect.with(stroke: none, width: 100%, height: 100%)
  #let gray_cell = rect.with(
    stroke: none,
    width: 110%,
    height: 100%,
    fill: gray,
    radius: (rest: 10pt),
  )
  #grid(
    columns: (auto, 52%),
    column-gutter: 5pt,
    cell[- DMI理论

      我们将DMI能量离散化
      $ E_("DM") = bold(D_(i,j)) dot sum_(i,j)(bold(m_i) times bold(m_j)) $
      通过场和能量的关系$H=- nabla_m E$可以得到离散化的DMI有效场的表达式
      $ attach(limits(H), t: arrow, br: "DM") = - 1/mu_0 (diff E)/(diff arrow(m)_i)= -D/(mu_0 M_s)sum(attach(limits(m), t: arrow)_j times attach(limits(r), t: arrow)_(i,j)) $

    ],
    gray_cell[
      #set text(size: 13pt)
      #show raw.where(block: true): set par(justify: false)
      ```python
     extern void calcai_3DMIfrac(int *nmx, int *nmy, int *nmz,
                                     unsigned char type[],
                                     double D_nu[256],
                                     double D_mix_nu[256],
                                     double Dx[], double Dy[],
                                     double Dz[], double D_dir[3],
                                     double Dc_dir[3]);
                                     //计算面内以及层间DMI系数设置DMI方向
         extern double ex_dmi(int *nmx, int *nmy, int *nmz, double *dx,
                             double *dy,double *dz, double Dx[],
                             double Dy[], double Dz[],
                              double mn[], double m[][3]);
                               //计算体系DMI能量
         void calc_Hd2(double m[][3], double mn[], double Hx[][3],
                       double heff[][3],
                       int ipmax, int *nx0, int *ny0, int *nz0, int idx,
                       int idy,
                       int idz, int iaxdx, int iaxdy, int iaxdz, int iaydx,
                       int iaydy,
                       int iaydz, int iazdx, int iazdy, int iazdz,
                       double *dx,
                       double *dy, double *dz, double *Dx, double *Dy,
                       double *Dz);//计算DMI有效场

                                                      ```

    ],
  )

]

#slide(
  title: "自旋轨道矩诱导的反铁磁磁矩翻转",
  new-section: "退磁场计算并行优化",
)[
  #let cell = rect.with(stroke: none, width: 94%, height: 100%)
  #let gray_cell = rect.with(
    stroke: none,
    width: 110%,
    height: 100%,
    fill: gray,
    radius: (rest: 10pt),
  )
  #set text(size: 18pt)
  #grid(
    columns: 2,
    column-gutter: 5pt,
    cell[
      #set align(left + horizon)
      在通过静磁能求解退磁能时，其离散形式为
      $ E_d = - mu_0 /2 sum_(i,j,k=1,2,dots)H_(d,i,j,k) dot M_(i,j,k) delta V $
      我们通过求解Poisson方程的方法来获得磁标势$V$，$(nabla^2 V_d)_(i,j,k)=- rho_(i,j,k)$。采用松弛迭代法(SOR)求解Poisson方程，SOR方法主要是通过求解
      $ underbrace(A x=b, "linear equation") =>(D+ omega L)x = omega b -[omega U+(omega -1)D]x $其中D为A的对角矩阵，L为上角矩阵，U为下角矩阵，$omega$为弛豫因子。
    ],
    gray_cell[
      矩阵非常适合用CUDA计算，可以降低时间和空间复杂度，提高计算速度。
      #text(
        size: 15pt,
        ```c
__global__ void updateVKernel_even(double *d_V, double *d_chrg, int *d_maxiter,int ndx, int ndy, int ndz, int idx, int idy, int idz, int icdx, int icdy, int icdz, double dx2, double dy2, double dz2, double diag, double w,double dVn,double dV,double resid,double maxf,double minf,double medf,double medf2);{
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= 1 && ix < ndx && iy >= 1 && iy < ndy && iz >= 1 && iz < ndz && (ix + iy + iz) % 2 == 1){
       ...
    }
}
```
      )
      可以将原本三层嵌套的for循环，通过对线程的索引降低为一次并行计算。
    ],
  )
]

#slide(
  title: "自旋轨道矩诱导的反铁磁磁矩翻转",
  new-section: "程序验证",
)[
  #set text(size: 20pt)
  #let gray_cell = rect.with(
    stroke: none,
    width: 100%,
    height: 80%,
    fill: gray,
    radius: (rest: 10pt),
  )
  #let cell = rect.with(stroke: none, width: 100%, height: 80%)
  - 我们利用层间DMI作用实现磁矩无场翻转作为程序验证的案例。
  #grid(
    columns: (35%, auto, 28%),
    gray_cell[
      - 程序参数
      #set align(left)

      $ M_s &= 1.5 times 10^6 A\/m \
      D_("inter") &= 2.0 times 10^(-5)J\/m^2 \
      J_("inter") &= 2.8 times 10^(-5) J\/m^2 \
      K &= 1.5 times 10^6 J\/m^3 \
      "cell size" &= 1 times 1 times 4 "nm" \
      "cell num" &=2 times 50 times 50 \ $
    ],
    cell[
      #set text(size: 17pt)
      #figure(image("pic/sot_switching.png"), caption: [层间DMI实现无场翻转结构])<glacier>

      @glacier 展示了$#ce("Pt/Co/Pt/Co/Pt")$结构的示意图，两个Co层之间有RKKY和层间DMI效应。顶部和底部$#ce("Pt")$层充当导体，用于施加电流pulse并生成
      SOT。

    ],
    cell[- 物理过程
      #set text(size: 18pt)
      #fletcher.diagram(
        cell-size: 10pt,
        spacing: 2.5em,
        node-stroke: black,
        node-fill: blue.lighten(90%),
        node((0, 3), "施加电流pulse产生SOT"),
        node((0, 2), "磁矩发生翻转"),
        node((0, 1), "撤去电流使磁矩弛豫"),
        node((0, 0), [在DMI作用下磁矩发生翻转 \ 翻转方向取决去DM常数方向]),
        edge((0, 0), (0, 1), "<="),
        edge((0, 1), (0, 2), "<="),
        edge((0, 2), (0, 3), "<="),
        // edge((1,0), (0,0), "..>", bend: -0deg),
      )
    ],
  )
]

#slide(title: "自旋轨道矩诱导的反铁磁磁矩翻转",new-section:"验证数据")[
#let cell = rect.with(stroke: none, width: 100%, height: 100%)
#let graycell = rect.with(stroke: none, fill: gray, radius: (rest: 10pt),height:90%,width:100%) 
#grid(columns:2,column-gutter:5pt,
cell[
  #figure(image("pic/data_unfield_switching.png",height:80%), caption: [层间DMI实现无场翻转物理过程])<figure3>
  ],
graycell[
  #set align(horizon)
  #set text(size: 18pt)
   - @figure3 (`b`)展示了在前10ns未加电流时，FM结构中磁矩的初始状态。@figure3 (`c`)展示了在10ns前AFM结构中磁矩的初始状态。

   - 随着10ns后施加电流pulse，在SOT效应的作用下，磁矩发生翻转，倒向了面内方向。

   - 电流pulse结束后，在DMI作用下，系统弛豫足够长时间后，磁矩在层间DMI的作用下发生了可控翻转，翻转方向取决于DMI方向。
   
    ])
  // ],
// )
]

// #slide(
//   title: "研究成果",
//   new-section: "论文情况",
// )[
//   #set text(size: 18pt)
//   #set align(horizon)
//   #let graycell = rect.with(stroke: none, fill: gray, radius: (rest: 10pt), height: 70%)
//   #v(2em)
//   #graycell[
//     - *目前有一篇关于层间DMI实现无场翻转的二作论文在投*
//
//     [1]Zheng, Cuixiu, Wenqing He, Mo Zhu, Caihua Wan, Xiufeng Han和Yaowen Liu.
//     《Interlayer Dzyaloshinskii–Moriya Interaction Assisted Field-Free Magnetization
//     Switching》
//   ]
// ]

// #slide(title: "下一步研究方向", new-section: "DMI对激发Thz频率的影响")[
//   #set text(size: 18pt)
//   #let graycell= rect.with(stroke: none, fill: gray, radius: (rest: 10pt))
//   #let cell= rect.with(stroke:none,height:55%)
//   #let fullcell= rect.with(stroke:none,width:110%)
// // == 我们使用`Vampire`作为我们原子尺度反铁磁模拟软件
// == *在反铁磁体系中探究DMI对电流诱导自旋波激发过程的影响*
// #v(1em)
// - 写出系统的哈密顿量
// #grid(rows: 2, 
// graycell[
//   $ H_s = J_1 sum_(<i,j>_("xy")) bold(S)_i dot bold(S)_j + J_2 sum_(<i,j>_z) bold(S)_i dot bold(S)_j +overbrace(sum_(<i,j>_("xy")) D_(i,j) dot bold(S)_i times bold(S)_j,D_("ij")^k=D^k dot (bold(z) times bold( u_("ij") ))) - sum_i K (accent(n,hat)_i dot bold(S)_i)^2 $],
// cell[
//   #grid(columns:2,
//   fullcell[
//   #figure(image("pic/layer_DMI.png",height:80%), caption: [`R A Duine, et al., Nature Physics, 14(2018), 217-219.`
//  \ 层间DMI示意图])
//   ],
//   fullcell[
//     #v(3em)
//     由于我们的体系是重金属层和反铁磁层双层结构，由于界面处的对称性破缺，磁性原子对和非磁性原子之间会产生较大的DMI。
//   ]
// )]


// )


// ]

#focus-slide[
  #set align(center + horizon)
  欢迎各位老师批评指正！

]



// #matrix-slide[
//   left
// ][
//   middle
// ][
//   right
// ]

// #matrix-slide(columns: 1)[
//   top
// ][
//   bottom
// ]:SymbolsOutline<CR>

// #matrix-slide(columns: (1fr, 2fr, 1fr), ..(lorem(8),) * 9)
