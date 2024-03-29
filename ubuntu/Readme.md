## bash command
### Useful commands
- `adduser [username]`
- ` ps -aux `
- ` df -h `
- ` du -sh ./* `

### Useful packages
- [CPU Temperature Monitor](https://askubuntu.com/questions/15832/how-do-i-get-the-cpu-temperature)
### Start script and service when boot
- [부팅 시 자동 시작되는 스크립트](https://nobilitycat.tistory.com/entry/%EB%A6%AC%EB%88%85%EC%8A%A4-%EC%8B%9C%EC%9E%91-%EC%8B%9C-%EC%9E%90%EB%8F%99%EC%9C%BC%EB%A1%9C-%EC%8B%A4%ED%96%89-%EB%90%A0-%ED%94%84%EB%A1%9C%EA%B7%B8%EB%9E%A8-%EB%93%B1%EB%A1%9D%ED%95%98%EA%B8%B0)
```bash
cd /etc/init.d
touch start.sh
cat "nohup jupyter lab &" >> start.sh
chmod 755 start.sh
update-rc.d start.sh defaults
```

## 하드디스크 추가
 - [우분투 하드 추가](http://blog.naver.com/PostView.nhn?blogId=jonghong0316&logNo=140146721316)
 - [우분투 적출 하드 추가 + 권한설정](https://juheo.tistory.com/7)

## [R 최신버전 설치](https://wikidocs.net/52543)


## Pyenv

## 도커
- [도커 설치](https://blog.cosmosfarm.com/archives/248/%EC%9A%B0%EB%B6%84%ED%88%AC-18-04-%EB%8F%84%EC%BB%A4-docker-%EC%84%A4%EC%B9%98-%EB%B0%A9%EB%B2%95/)
- [우분투 도커 설치 using 캐글 이미지 with CPU/GPU 버전](https://teddylee777.github.io/linux/docker%EB%A5%BC-%ED%99%9C%EC%9A%A9%ED%95%98%EC%97%AC-%EB%94%A5%EB%9F%AC%EB%8B%9D-%ED%99%98%EA%B2%BD%EA%B5%AC%EC%84%B1.md)


## pyenv 설치
 - [pyenv 설치](https://nachwon.github.io/pyenv-virtualenv/)


## Jupyter

- `jupyter kernelspec list`
- `nohup jupyter notebook &`
- `jupyter-lab list`
- `jupyter --paths`
- `jupyter notebook stop [port number]`
- `lsof -i :8888`
- `kill -9 $(lsof -i :8888 | sed -n 2p | awk -F ' ' '{print $2}')`


### [Jupyter notebook에 R 추가하기](https://yahwang.github.io/posts/27)




## LightGBM with GPU
https://www.kaggle.com/kirankunapuli/ieee-fraud-lightgbm-with-gpu

### INSTALL NVIDIA GRAPHIC DRIVER
- https://inpages.tistory.com/149
### cmake에서 libOpenCL.so 대신 libOpenCL.so.1 사용하기
- cmake -DUSE_GPU=1 -DOpenCL_LIBRARY=/usr/local/cuda-11.1/lib64/libOpenCL.so.1 -DOpenCL_INCLUDE_DIR=/usr/local/cuda-11.1/include/ ..
- https://stackoverflow.com/questions/50627963/lightgbmerror-bno-opencl-device-found


```
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
```

