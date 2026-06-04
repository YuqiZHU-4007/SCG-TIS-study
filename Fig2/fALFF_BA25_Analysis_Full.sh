#!/bin/bash

# ==============================================================================
# BA25 脑区 fALFF 交互作用分析：全流程自动化脚本 (补全 ACF 估计步骤)
# ==============================================================================

# 1. 设置基础路径与环境变量
BASE_DIR="/mnt/e/fALFF_data/AFNI_fALFF"
MASK_FILE="${BASE_DIR}/BA25_L_on_residGrid.nii"
DATA_TABLE="${BASE_DIR}/LMEr_fALFF_noTime.txt"

# 切换到工作目录
cd ${BASE_DIR} || exit

echo "----------------------------------------------------------------"
echo "第一阶段：线性混合效应模型（LME）建模"
echo "----------------------------------------------------------------"
# 生成统计图及残差文件，残差文件用于后续 ACF 估计[cite: 1, 5]
3dLMEr -prefix LMEr_fALFF_GroupXPhase_ver.2 \
       -resid LMEr_fALFF_GroupXPhase_ver.2_resid \
       -jobs 8 \
       -model 'Group*Phase+(1|Subj)' \
       -gltCode Phase_stim_vs_pre Phase : '1*stim' '-1*pre' \
       -gltCode GroupXPhase Group : '1*10Hz' '-1*Sham' Phase : '1*stim' '-1*pre' \
       -dataTable @${DATA_TABLE}

echo "----------------------------------------------------------------"
echo "第二阶段：估计空间自相关参数 (ACF)"
echo "----------------------------------------------------------------"
# 使用 3dFWHMx 从残差文件中计算 ACF 参数
# -mask: 限制在 BA25 掩码内进行估计
# -acf: 计算 ACF 参数 (a, b, c)
# -input: 使用第一步生成的残差文件
3dFWHMx -mask ${MASK_FILE} \
        -acf \
        -input LMEr_fALFF_GroupXPhase_ver.2_resid+orig > BA25_acf_params.txt

# 从输出文件中提取最后一行的 3 个 ACF 数值
# 之前计算的值为：0.993985 3.56403 0.425826
ACF_VALS=$(tail -n 1 BA25_acf_params.txt)
echo "提取到的 ACF 参数为: ${ACF_VALS}"

echo "----------------------------------------------------------------"
echo "第三阶段：初步统计阈值化 (p < 0.005)"
echo "----------------------------------------------------------------"
# 提取交互作用 Z 分数并进行双端硬阈值处理[cite: 1, 5]
3dcalc -a "LMEr_fALFF_GroupXPhase_ver.2+orig[6]" \
       -expr 'a*step(abs(a)-1.96)' \
       -prefix GroupXPhase_ver.3_Z_p005

echo "----------------------------------------------------------------"
echo "第四阶段：基于 ACF 的聚类校正模拟 (3dClustSim)"
echo "----------------------------------------------------------------"
# 使用提取到的动态参数进行模拟
3dClustSim -mask ${MASK_FILE} \
           -acf ${ACF_VALS} \
           -pthr 0.005 \
           -iter 10000 \
           -prefix BA25_clustsim

echo "----------------------------------------------------------------"
echo "第五阶段：聚类筛选与结果提取"
echo "----------------------------------------------------------------"
# 筛选满足簇大小（13 体素）的显著区域并生成报告[cite: 1, 4]
3dClusterize -inset GroupXPhase_ver.3_Z_p005+orig \
             -ithr 0 \
             -bisided -1.96 1.96 \
             -NN 2 \
             -mask ${MASK_FILE} \
             -clust_nvox 13 \
             -pref_map GroupXPhase_ver.3_Z_p005_mask2 \
             > GroupXPhase_mask2_Z_ver.3.txt

# 转换为 NIfTI 格式
3dAFNItoNIFTI -prefix GroupXPhase_ver.3_Z_p005_BA25.nii.gz GroupXPhase_ver.3_Z_p005_mask2+orig

echo "================================================================"
echo "分析流程结束。"
echo "================================================================"