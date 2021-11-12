/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#ifndef _CLASS2_H
#define _CLASS2_H

#include<iostream>
#include<cstring>	//memset�Cmemcpy
#include<new>		//��O����

class packer{
private:
	unsigned int max_size;	//�m�ۂ���o�C�g���Copteron�̏ꍇint=long int(32bit=2G)�����C�܂�������Dlongna�Ȃ�64bitOS���������͂Ȃ������...9223372036854775808
	int* index;		//bin�̓�4�o�C�g��index�Ɏg��	//�p�b�L���O���Ă���o�C�g���̕ێ�
	char* bin;		//�{�́Cchar�Ȃ͎̂���ȁHunsigned��t�������??
public:
	packer(){
		max_size =sizeof(int);
		try{
			bin = new char[max_size];
		}catch(const std::bad_alloc&){
			std::cerr<<"allocerr"<<std::endl;
			exit(1);
		}
		memset(bin,0,max_size*sizeof(char));
		index=(int*)bin;
		*index=sizeof(int);
	}
	packer(const unsigned int n){//�m�ۃo�C�g�����
		max_size = n;
		if(n<sizeof(int))max_size =sizeof(int);
		try{
			bin = new char[max_size];
		}catch(const std::bad_alloc&){
			std::cerr<<"allocerr"<<std::endl;
			exit(1);
		}
		memset(bin,0,max_size*sizeof(char));
		index=(int*)bin;
		*index=sizeof(int);
	}
	~packer(){
		index=NULL;
		max_size =0;
		delete[] bin;
	}
	
	
 	char* operator[](const unsigned int i){return bin+i;}//void*�ɂ��ׂ���//����킵���ȁCi�Ԗڂ̗v�f�A�N�Z�X�łȂ��Ci�o�C�g�ڂ̗v�f�A�N�Z�X
	
	template <class T>void pack(const T* const in){
		if(*index+sizeof(T)<max_size){
			memcpy(bin+*index,in,sizeof(T));
			*index+=sizeof(T);
		}else{
			resize(*index+sizeof(T)+1024);
			pack(in);//�ċA
		}
	}
	template <class T>void pack(const T* const in,const int n){//"�A��������"�z��	�A������Ȃ��Ƒ�ςȂ��ƂɂȂ�܂���
		if(*index+n*sizeof(T)<max_size){
			memcpy(bin+*index,in,n*sizeof(T));
			*index+=n*sizeof(T);
		}else{
			resize(*index+n*sizeof(T)+1024);
			pack(in,n);
		}
	}
	
	void resize(const unsigned int i){
		if(max_size<i){//�k���͂����܂���D
			char* tmp;
			try{
				tmp=new char[i];
			}catch(const std::bad_alloc&){
				std::cout<<"allocerr"<<std::endl;
				exit(1);
			}
			memset(tmp, 0 , i * sizeof(char));
			memcpy(tmp,bin,max_size*sizeof(char));
			delete[] bin;
			bin=tmp;
			max_size=i;
			index=(int*)bin;
		}
	}
	void clear(){
		memset(bin,0,*index*sizeof(char));
		*index=sizeof(int);
	}
	void clear(const int bin_size){
		memset(bin,0,bin_size*sizeof(char));
		*index=sizeof(int);
	}
	void allclear(){
		memset(bin,0,max_size*sizeof(char));
		*index=sizeof(int);
	}
	int size()const{return *index;}
	unsigned int maxsize()const{return max_size;}
//	void size(const int n){index=n;}
	
};




#endif


#if 0
class packer
��MPI�Ńf�[�^�𑗂鎞�Ƀf�[�^���p�b�P�[�W���O���邽�߂̃N���X�ł��D
��{�\���͓�����Ƃ͂��Ă��Ȃ��̂œǂ߂΂킩��͂��D
unpack�֐�������Ă��Ȃ��̂͒P�Ȃ�ӑĂł��D
pack-unpack�����܂������������
�����Ǝg���₷���͂Ȃ肻�����ȁD���Ȃ����ǁD
MPI�̗��_�����\�������Ă���̂ŁC���܂��ł�����@���ʂɂ��邩������Ȃ��ł��D
�܂����ł�c++�̗��_�������薳�����Ă��܂����D

�g�����Ƃ��ẮC
���Ɍ^���ʗp�̎��ʎq��u���Ă��̎���

#endif//�R�����g

