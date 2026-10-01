#!/bin/bash
# ABodyBuilder2 Docker 运行脚本
# ============================================================
# 确保 Docker 已安装，NVIDIA Container Toolkit 已配置

# 创建数据目录
mkdir -p /tmp/abodybuilder2_data

# 写入 FASTA 文件
cat > /tmp/abodybuilder2_data/VH_VL.fasta << 'FASTA'
>IA6-13-8_VH
EVQLQQSGPELVKPGTSVKVSckasGYSFTDYNiYWVKQSHGKSLEWIGIDPYNDDTSYNQKFKGKATLTVDKSSSTAFMHLNSLTSEDSAVYYCTRWEITDPWYFDVWGAGTTVTVSS
>IA6-13-8_VL
DIVMTQSHKFMSISVGDRASIMckasQDVGTAvsWYQQKPGQSPKALTYWASRRHTGVPDRFTGSGSGTDFTLIIGNVQEDLAAYFCQQYSSYPYTFGGGTKLEIK
FASTA

# 运行 Docker 容器
docker run -it \
  -v /tmp/abodybuilder2_data:/data \
  oxpig/abodybuilder2:latest \
  ABodyBuilder2 \
  --fasta_file /data/VH_VL.fasta \
  --output_dir /data/output

# 复制结果到当前目录
cp -r /tmp/abodybuilder2_data/output/* /home/xinlab/home/ws/phase2/results/11_antibody_structure/

echo "✅ 结果已保存到: /home/xinlab/home/ws/phase2/results/11_antibody_structure"
