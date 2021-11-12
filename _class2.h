/*
//Brownian Dynamics simulation of filaments
//Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
*/
#ifndef _CLASS2_H
#define _CLASS2_H

#include<iostream>
#include<cstring>	//memset，memcpy
#include<new>		//例外処理

class packer{
private:
	unsigned int max_size;	//確保するバイト数，opteronの場合int=long int(32bit=2G)だし，まぁいいや．longnaなら64bitOSだから上限はないけれど...9223372036854775808
	int* index;		//binの頭4バイトをindexに使う	//パッキングしているバイト数の保持
	char* bin;		//本体，charなのは趣味かな？unsignedを付けろって??
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
	packer(const unsigned int n){//確保バイト数代入
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
	
	
 	char* operator[](const unsigned int i){return bin+i;}//void*にすべきか//紛らわしいな，i番目の要素アクセスでなく，iバイト目の要素アクセス
	
	template <class T>void pack(const T* const in){
		if(*index+sizeof(T)<max_size){
			memcpy(bin+*index,in,sizeof(T));
			*index+=sizeof(T);
		}else{
			resize(*index+sizeof(T)+1024);
			pack(in);//再帰
		}
	}
	template <class T>void pack(const T* const in,const int n){//"連続メモリ"配列	連続じゃないと大変なことになりますね
		if(*index+n*sizeof(T)<max_size){
			memcpy(bin+*index,in,n*sizeof(T));
			*index+=n*sizeof(T);
		}else{
			resize(*index+n*sizeof(T)+1024);
			pack(in,n);
		}
	}
	
	void resize(const unsigned int i){
		if(max_size<i){//縮小はさせません．
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
はMPIでデータを送る時にデータをパッケージングするためのクラスです．
基本構造は難しいことはしていないので読めばわかるはず．
unpack関数を作っていないのは単なる怠惰です．
pack-unpackをうまく書き換えれば
もっと使いやすくはなりそうかな．やらないけど．
MPIの利点を結構無視しているので，うまくできる方法が別にあるかもしれないです．
まぁついでにc++の利点もがっつり無視していますが．

使い方としては，
頭に型判別用の識別子を置いてその次に

#endif//コメント

