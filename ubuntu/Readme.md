# Pyenv

# 도커

- [우분투 도커 설치 using 캐글 이미지 with CPU/GPU 버전](https://teddylee777.github.io/linux/docker%EB%A5%BC-%ED%99%9C%EC%9A%A9%ED%95%98%EC%97%AC-%EB%94%A5%EB%9F%AC%EB%8B%9D-%ED%99%98%EA%B2%BD%EA%B5%AC%EC%84%B1.md)

# Jupyter

# LightGBM with GPU
https://www.kaggle.com/kirankunapuli/ieee-fraud-lightgbm-with-gpu

## cmake에서 libOpenCL.so 대신 libOpenCL.so.1 사용하기
- cmake -DUSE_GPU=1 -DOpenCL_LIBRARY=/usr/local/cuda-11.1/lib64/libOpenCL.so.1 -DOpenCL_INCLUDE_DIR=/usr/local/cuda-11.1/include/ ..
- https://stackoverflow.com/questions/50627963/lightgbmerror-bno-opencl-device-found


₩₩₩
import lightgbm as lgb
import seaborn as sns

iris = sns.load_dataset('iris')
iris.head()

%%bash
cd /home/png/Install/
cd LightGBM
rm -r build
mkdir build
cd build
cmake -DUSE_GPU=1 -DOpenCL_LIBRARY=/usr/local/cuda-11.1/lib64/libOpenCL.so.1 -DOpenCL_INCLUDE_DIR=/usr/local/cuda-11.1/include/ ..
make -j$(nproc)

%%bash
python3 -m pip install wheel setuptools


sudo -S python3 /home/png/Install/LightGBM/python-package/setup.py install --precompile < /home/png/passwd/png

sudo -S mkdir -p /etc/OpenCL/vendors
echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd


d_train= lgb.Dataset(iris.loc[:,iris.columns != "species"], iris.species)

params = {
        'objective': 'multiclass',
        'device': 'gpu', 
        # 'num_boost_round': n_estimators,
        "metric": "multi_logloss", "num_class": 3
    }
    
    
model = lgb.train(
                        params,
                        train_set       = d_train,
                        num_boost_round = 100,
                        valid_sets =[d_train],
                        valid_names=['train'],
                        verbose_eval    = 200, 
                        early_stopping_rounds = 10
                        )
₩₩₩     
