#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define mod(a) ((a)<0) ? (-a) : (a)
typedef struct
{
    unsigned char B,G,R;
} RGB;

typedef struct
{
    int X1,Y1,X2,Y2;
    int cifra_det;
    double valoare_corelatie;
} fereastra;

typedef struct
{
    fereastra *V;
    unsigned int NumarDetectii;
} ferestre_sablon;

unsigned int xorshift32(unsigned int *seed)
{
    unsigned int x;
    x=*seed;
    x^= x << 13;
    x^= x >> 17;
    x^= x << 5;
    *seed=x;
    return x;
}

void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie)
{
    FILE *fin, *fout;
    unsigned int dim_img, latime_img, inaltime_img;
    unsigned char *pRGB, aux;
    pRGB=(unsigned char*)malloc(3*sizeof(unsigned char));
    if(pRGB==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia grayscale_image./n");
        return ;
    }
    fin = fopen(nume_fisier_sursa, "rb");
    if(fin == NULL)
    {
        printf("Nu am gasit imaginea sursa din care citesc pentru a transforma imaginea in grayscale.\n");
        return;
    }

    fout = fopen(nume_fisier_destinatie, "wb+");

    fseek(fin, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fin);

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

    //copiaza octet cu octet imaginea initiala in cea noua
    fseek(fin,0,SEEK_SET);
    unsigned char c;
    while(fread(&c,1,1,fin)==1)
    {
        fwrite(&c,1,1,fout);
        fflush(fout);
    }
    fclose(fin);

    //calculam padding-ul pentru o linie
    int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    fseek(fout, 54, SEEK_SET);
    int i,j;
    for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++)
        {
            //citesc culorile pixelului
            fread(pRGB, 3, 1, fout);
            //fac conversia in pixel gri
            aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
            pRGB[0] = pRGB[1] = pRGB[2] = aux;
            fseek(fout, -3, SEEK_CUR);
            fwrite(pRGB, 3, 1, fout);
            fflush(fout);
        }
        fseek(fout,padding,SEEK_CUR);
    }
    fclose(fout);
    free(pRGB);
}

RGB* IncarcareFormaLiniarizata(char *CaleImagine)
{
    FILE *f;
    unsigned int Width,Height,padding;
    RGB *V;
    int i,j;
    f=fopen(CaleImagine,"rb");
    if(f==NULL)
    {
        printf("Eroare la deschiderea imaginii pentru a o liniariza.\n");
        return NULL;
    }
    fseek(f,18,SEEK_SET);
    fread(&Width,sizeof(unsigned int),1,f);
    fread(&Height,sizeof(unsigned int),1,f);
    V=(RGB*)malloc(Width*Height*sizeof(RGB));
    if(V==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia IncarcareFormaLiniarizata.\n");
        return NULL;
    }
    if(Width%4)
        padding=4-(Width*3)%4;
    else
        padding=0;
    fseek(f,54,SEEK_SET);
    for(i=Height-1; i>=0; i--)
    {
        for(j=0; j<Width; j++)
            fread((V+Width*i+j),3,1,f);
        fseek(f,padding,SEEK_CUR);
    }
    fclose(f);
    return V;
}

unsigned int* generare_permutare_random(unsigned int *R,unsigned int nr)
{
    unsigned int i,r,*permutare,aux;
    permutare=(unsigned int *)malloc(nr*sizeof(unsigned int));
    if(permutare==NULL)
    {
        printf("Eroare la alocarea dinamica a memoriei in functia generare_permutare_random.\n");
        return NULL;
    }
    for(i=0; i<nr; i++)
        *(permutare+i)=i;
    for(i=nr-1; i>=1; i--)
    {
        r=*(R+nr-i-1);
        if(r>i)
            r=r%(i+1);
        aux=*(permutare+r);
        *(permutare+r)=*(permutare+i);
        *(permutare+i)=aux;
    }
    return permutare;
}

void SalvareExternFormaLiniarizata(char *cale_imagine,char *cale_imagine_destinatie,RGB* VectorImagine)
{
    FILE *f,*g;
    int i,j;
    unsigned int Width,Height,padding;
    unsigned char *header,val=0;
    header=(unsigned char*)malloc(54);
    if(header==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia SalvareExternFormaLiniarizata.\n");
        return ;
    }
    f=fopen(cale_imagine,"rb");
    if(f==NULL)
    {
        printf("Eroare la deschiderea imaginii sursa.\n");
        return;
    }
    fread(header,54,1,f);
    fseek(f,18,SEEK_SET);
    fread(&Width,sizeof(unsigned int),1,f);
    fread(&Height,sizeof(unsigned int),1,f);
    fclose(f);
    if(Width%4)
        padding=4-(Width*3)%4;
    else
        padding=0;
    g=fopen(cale_imagine_destinatie,"wb");
    if(g==NULL)
    {
        printf("Eroare la deschiderea fisierului destinatie a imaginii.\n");
        return;
    }
    fwrite(header,54,1,g);
    for(i=Height-1; i>=0; i--)
    {
        fwrite((VectorImagine+i*Width),3,Width,g);
        if(padding)
            for(j=1; j<=padding; j++)
                fwrite(&val,1,1,g);
    }
    fclose(g);
    free(header);
}

void criptare(char *cale_imagine,char *cale_imagine_criptata,char *cheia_secreta)
{
    FILE *f,*g;
    RGB *P,*P_permutat;
    unsigned int *R,Width,Height,seed,SW,*Permutare,i;
    f=fopen(cale_imagine,"rb");
    if(f==NULL)
    {
        printf("Eroare la deschiderea imaginii sursa.\n");
        return;
    }
    fseek(f,18,SEEK_SET);
    fread(&Width,sizeof(unsigned int),1,f);
    fread(&Height,sizeof(unsigned int),1,f);
    fclose(f);
    g=fopen(cheia_secreta,"r");
    if(g==NULL)
    {
        printf("Eroare la deschiderea fisierului cu cheia secreta.\n");
        return;
    }
    fscanf(g,"%u",&seed);
    fscanf(g,"%u",&SW);
    fclose(g);
    P=IncarcareFormaLiniarizata(cale_imagine); //am liniarizat imaginea initiala
    if(P==NULL)
    {
        printf("Eroare la liniarizarea imaginii.\n");
        return;
    }

//generare 2*Width*Height-1 numere aleatoare pornind de la valoarea seed
    R=(unsigned int*)malloc((2*Width*Height-1)*sizeof(unsigned int));
    if(R==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia de criptare.\n");
        return ;
    }
    for(i=0; i<2*Width*Height-1; i++)
        *(R+i)=xorshift32(&seed);

//generare permutare random cu ajutor primilor Width*Height-1 numere generate mai sus
    Permutare=generare_permutare_random(R,Width*Height);
    if(Permutare==NULL)
        return;
//permutam pixelii imaginii liniarizate conform permutarii generate mai sus
    P_permutat=(RGB*) malloc((Width*Height)*sizeof(RGB));
    if(P_permutat==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia de criptare.\n");
        return ;
    }
    for(i=0; i<Width*Height; i++)
        *(P_permutat + *(Permutare+i)) = *(P+i);

//XOR-am pixelii imaginii liniarizate
    unsigned char *aux1,*aux2;
    aux1=(unsigned char*)&SW;
    aux2=(unsigned char*)(R+Width*Height-1);
    P_permutat->R = (*(aux1+2)) ^ (P_permutat->R) ^ (*(aux2+2));
    P_permutat->G = (*(aux1+1)) ^ (P_permutat->G) ^ (*(aux2+1));
    P_permutat->B = (*aux1) ^ (P_permutat->B) ^ (*aux2);
    for(i=0; i<Width*Height-1; i++)
    {
        aux2=(unsigned char*)(R+Width*Height+i);
        (P_permutat+i+1)->R = ((P_permutat+i)->R) ^ ((P_permutat+i+1)->R) ^ (*(aux2+2));
        (P_permutat+i+1)->G = ((P_permutat+i)->G) ^ ((P_permutat+i+1)->G) ^ (*(aux2+1));
        (P_permutat+i+1)->B = ((P_permutat+i)->B) ^ ((P_permutat+i+1)->B) ^ (*aux2);
    }
    SalvareExternFormaLiniarizata(cale_imagine,cale_imagine_criptata,P_permutat);
    free(R);
    free(P_permutat);
    free(Permutare);
    free(P);
}

void decriptare(char *cale_imagine_criptata,char *cale_imagine_decriptata,char *cheia_secreta)
{
    FILE *f,*g;
    RGB *P,*P_permutat;
    unsigned int *R,Width,Height,seed,SW,*Permutare,i;
    unsigned char *aux1,*aux2;

//citirea din fisier
    f=fopen(cale_imagine_criptata,"rb");
    if(f==NULL)
    {
        printf("Eroare la deschiderea imaginii sursa.\n");
        return;
    }
    fseek(f,18,SEEK_SET);
    fread(&Width,sizeof(unsigned int),1,f);
    fread(&Height,sizeof(unsigned int),1,f);
    fclose(f);
    g=fopen(cheia_secreta,"r");
    if(g==NULL)
    {
        printf("Eroare la deschiderea fisierului cu cheia secreta.\n");
        return;
    }
    fscanf(g,"%u",&seed);
    fscanf(g,"%u",&SW);
    fclose(g);

//generare 2*Width*Height-1 numere aleatoare pornind de la valoarea seed
    R=(unsigned int*)malloc((2*Width*Height-1)*sizeof(unsigned int));
    if(R==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia de decriptare.\n");
        return ;
    }
    for(i=0; i<2*Width*Height-1; i++)
        *(R+i)=xorshift32(&seed);

    P=IncarcareFormaLiniarizata(cale_imagine_criptata);
    if(P==NULL)
    {
        printf("Eroare la liniarizarea imaginii.\n");
        return;
    }
//XOR-am pixelii imaginii liniarizate pentru a ajunge la valorile initiale ale pixelilor

    for(i=Width*Height-1; i>=1; i--)
    {
        aux2=(unsigned char*)(R+Width*Height+i-1);
        (P+i)->R = ((P+i-1)->R) ^ ((P+i)->R) ^ (*(aux2+2));
        (P+i)->G = ((P+i-1)->G) ^ ((P+i)->G) ^ (*(aux2+1));
        (P+i)->B = ((P+i-1)->B) ^ ((P+i)->B) ^ (*aux2);
    }
    aux2=(unsigned char*)(R+Width*Height-1);
    aux1=(unsigned char*)&SW;
    P->R = (*(aux1+2)) ^ (P->R) ^ (*(aux2+2));
    P->G = (*(aux1+1)) ^ (P->G) ^ (*(aux2+1));
    P->B = (*aux1) ^ (P->B) ^ (*aux2);

//generare permutare random cu ajutor primilor Width*Height-1 numere generate mai sus
    Permutare=generare_permutare_random(R,Width*Height);
    if(Permutare==NULL)
        return ;
//permutam pixelii imaginii liniarizate conform permutarii inverse
    P_permutat=malloc((Width*Height)*sizeof(RGB));
    if(P_permutat==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia de criptare.\n");
        return ;
    }
    for(i=0; i<Width*Height; i++)
        *(P_permutat + i) = *(P+*(Permutare+i));

    SalvareExternFormaLiniarizata(cale_imagine_criptata,cale_imagine_decriptata,P_permutat);
    free(R);
    free(P_permutat);
    free(P);
    free(Permutare);
}

void test_chi_patrat(char *cale_imagine)
{
    FILE *f;
    unsigned int Width,Height,i;
    f=fopen(cale_imagine,"rb");
    if(f==NULL)
    {
        printf("Eroare la deschiderea imaginii.\n");
        return;
    }
    fseek(f,18,SEEK_SET);
    fread(&Width,sizeof(unsigned int),1,f);
    fread(&Height,sizeof(unsigned int),1,f);
    fclose(f);
    RGB *Vector_imagine;
    unsigned int *Vector_frecvente_R,*Vector_frecvente_G,*Vector_frecvente_B;
    Vector_frecvente_R=(unsigned int*)calloc(256,sizeof(unsigned int));
    if(Vector_frecvente_R==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia test_chi_patrat.\n");
        return ;
    }
    Vector_frecvente_G=(unsigned int*)calloc(256,sizeof(unsigned int));
    if(Vector_frecvente_G==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia test_chi_patrat.\n");
        return ;
    }
    Vector_frecvente_B=(unsigned int*)calloc(256,sizeof(unsigned int));
    if(Vector_frecvente_B==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia test_chi_patrat.\n");
        return ;
    }
    double fmed,chi_patrat_R=0,chi_patrat_G=0,chi_patrat_B=0;
    fmed=(Width*Height)/256;
    Vector_imagine=IncarcareFormaLiniarizata(cale_imagine);
    if(Vector_imagine==NULL)
    {
        printf("Eroare la liniarizarea imaginii.\n");
        return;
    }
    for(i=0; i<Width*Height; i++)
    {
        *(Vector_frecvente_R+((Vector_imagine+i)->R))=*(Vector_frecvente_R+((Vector_imagine+i)->R))+1;
        *(Vector_frecvente_G+((Vector_imagine+i)->G))=*(Vector_frecvente_G+((Vector_imagine+i)->G))+1;
        *(Vector_frecvente_B+((Vector_imagine+i)->B))=*(Vector_frecvente_B+((Vector_imagine+i)->B))+1;
    }
    for(i=0; i<=255; i++)
    {
        chi_patrat_R+=((*(Vector_frecvente_R+i)-fmed)*(*(Vector_frecvente_R+i)-fmed))/fmed;
        chi_patrat_G+=((*(Vector_frecvente_G+i)-fmed)*(*(Vector_frecvente_G+i)-fmed))/fmed;
        chi_patrat_B+=((*(Vector_frecvente_B+i)-fmed)*(*(Vector_frecvente_B+i)-fmed))/fmed;
    }
    printf("(%lf,%lf,%lf)\n",chi_patrat_R,chi_patrat_G,chi_patrat_B);
    free(Vector_imagine);
    free(Vector_frecvente_B);
    free(Vector_frecvente_G);
    free(Vector_frecvente_R);
}

double mediavalorilor(RGB *V,unsigned int Width,unsigned int Height)
{
    int i;
    unsigned int n;
    double Smed=0;
    n=Width*Height;
    for(i=0; i<Height*Width; i++)
        Smed+=((V+i)->R);
    Smed=Smed/n;
    return Smed;
}

double deviatia(RGB *V,unsigned int Width,unsigned int Height,double med)
{
    int i;
    unsigned int n;
    double S=0;
    n=Width*Height;
    for(i=0; i<Height*Width; i++)
        S+=((V+i)->R-med)*((V+i)->R-med);
    S=S/(n-1);
    S=sqrt(S);
    return S;
}

double corelatie(RGB *S,unsigned int WidthS,unsigned int HeightS,RGB *I,unsigned int Width,unsigned int Height,int Iy,int Ix)
{
    RGB *ImagineSubSablon;
    ImagineSubSablon=(RGB *)malloc(WidthS*HeightS*sizeof(RGB));
    if(ImagineSubSablon==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia de calculare a corelatiei.\n");
        return -2;
    }
    int Sx,Sy,i,j,nr=0;
    Sx=Ix-(WidthS/2);
    Sy=Iy-(HeightS/2);

    for(i=Sy; i<Sy+(int)HeightS; i++)
        for(j=Sx; j<Sx+(int)WidthS; j++)
        {
            if(i>=0 && i<Height && j>=0 && j<Width)
                (ImagineSubSablon+nr)->R=(I+i*Width+j)->R;
            else
                (ImagineSubSablon+nr)->R=0;
            nr++;
        }
    double DS,DI,corel=0,Smed,Imed;
    Smed=mediavalorilor(S,WidthS,HeightS);
    Imed=mediavalorilor(ImagineSubSablon,WidthS,HeightS);
    DS=deviatia(S,WidthS,HeightS,Smed);
    DI=deviatia(ImagineSubSablon,WidthS,HeightS,Imed);
    for(i=0; i<HeightS; i++)
        for(j=0; j<WidthS; j++)
            corel+=(((ImagineSubSablon+i*WidthS+j)->R - Imed) * ((S+i*WidthS+j)->R - Smed)) / (DS * DI);
    corel=corel/(WidthS*HeightS);
    free(ImagineSubSablon);
    return corel;
}

ferestre_sablon template_matching(char *cale_imagine,char *cale_sablon,double ps)
{
    FILE *S,*I;
    RGB *VectorS,*VectorI;
    unsigned int WidthS,HeightS,Width,Height;
    int i,j,k;
    double rez;
    ferestre_sablon A;
    fereastra *aux;
    A.NumarDetectii=0;

    S=fopen(cale_sablon,"rb");
    if(S==NULL)
    {
        printf("Eroare la deschiderea imaginii sablon.\n");
        ferestre_sablon eroare;
        eroare.NumarDetectii=-1;
        return eroare;
    }
    fseek(S,18,SEEK_SET);
    fread(&WidthS,sizeof(unsigned int),1,S);
    fread(&HeightS,sizeof(unsigned int),1,S);
    fclose(S);

    I=fopen(cale_imagine,"rb");
    if(I==NULL)
    {
        printf("Eroare la deschiderea imaginii pe care se va opera functia de template matching.\n");
        ferestre_sablon eroare;
        eroare.NumarDetectii=-1;
        return eroare;
    }
    fseek(I,18,SEEK_SET);
    fread(&Width,sizeof(unsigned int),1,I);
    fread(&Height,sizeof(unsigned int),1,I);
    fclose(I);

    VectorS=IncarcareFormaLiniarizata(cale_sablon);
    if(VectorS==NULL)
    {
        printf("Eroare la liniarizarea imaginii.\n");
        ferestre_sablon eroare;
        eroare.NumarDetectii=-1;
        return eroare;
    }

    VectorI=IncarcareFormaLiniarizata(cale_imagine);
    if(VectorI==NULL)
    {
        printf("Eroare la liniarizarea imaginii.\n");
        ferestre_sablon eroare;
        eroare.NumarDetectii=-1;
        return eroare;
    }

    for(i=0; i<Height; i++)
        for(j=0; j<Width; j++)
        {
            rez=corelatie(VectorS,WidthS,HeightS,VectorI,Width,Height,i,j);
            if(rez>=ps)
            {
                aux=(fereastra *) malloc((A.NumarDetectii)*sizeof(fereastra));
                if(aux==NULL)
                {
                    printf("Eroare la alocarea memoriei dinamice in functia de template matching.\n");
                    ferestre_sablon eroare;
                    eroare.NumarDetectii=-1;
                    return eroare;
                }
                for(k=0; k<A.NumarDetectii; k++)
                {
                    (aux+k)->X1=(A.V+k)->X1;
                    (aux+k)->Y1=(A.V+k)->Y1;
                    (aux+k)->X2=(A.V+k)->X2;
                    (aux+k)->Y2=(A.V+k)->Y2;
                    (aux+k)->valoare_corelatie=(A.V+k)->valoare_corelatie;
                }
                free(A.V);
                A.NumarDetectii++;
                A.V=(fereastra *)malloc(A.NumarDetectii*sizeof(fereastra));
                if(A.V==NULL)
                {
                    printf("Eroare la alocarea memoriei dinamice in functia de template matching.\n");
                    ferestre_sablon eroare;
                    eroare.NumarDetectii=-1;
                    return eroare;
                }
                for(k=0; k<A.NumarDetectii-1; k++)
                {
                    (A.V+k)->X1=(aux+k)->X1;
                    (A.V+k)->Y1=(aux+k)->Y1;
                    (A.V+k)->X2=(aux+k)->X2;
                    (A.V+k)->Y2=(aux+k)->Y2;
                    (A.V+k)->valoare_corelatie=(aux+k)->valoare_corelatie;
                }
                (A.V+A.NumarDetectii-1)->X1=j-(WidthS/2);
                (A.V+A.NumarDetectii-1)->Y1=i-(HeightS/2);
                (A.V+A.NumarDetectii-1)->X2=j+(WidthS/2);
                (A.V+A.NumarDetectii-1)->Y2=i+(HeightS/2);
                (A.V+A.NumarDetectii-1)->valoare_corelatie=rez;
                free(aux);
            }
            else if(rez==-2) //va fi -2 daca nu se poate aloca suficienta memorie in functia de corelatie
            {
                ferestre_sablon eroare;
                eroare.NumarDetectii=-1;
                return eroare;
            }
        }
    free(VectorI);
    free(VectorS);
    return A;
}

void creare_tablou(ferestre_sablon *D,ferestre_sablon *Template,int cifra)
{
    int i;
    fereastra *aux;
    aux=(fereastra*) malloc(D->NumarDetectii * sizeof(fereastra));
    if(aux==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia de creare tablou cu toate detectiile.\n");
        return;
    }
    for(i=0; i<D->NumarDetectii; i++)
    {
        (aux+i)->valoare_corelatie=(D->V+i)->valoare_corelatie;
        (aux+i)->X1=(D->V+i)->X1;
        (aux+i)->Y1=(D->V+i)->Y1;
        (aux+i)->Y2=(D->V+i)->Y2;
        (aux+i)->X2=(D->V+i)->X2;
        (aux+i)->cifra_det=(D->V+i)->cifra_det;

    }
    D->NumarDetectii=(D->NumarDetectii) + (Template->NumarDetectii);
    free(D->V);
    D->V=(fereastra *)malloc(D->NumarDetectii * sizeof(fereastra));
    if(D->V==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice in functia de creare tablou cu toate detectiile.\n");
        return;
    }
    for(i=0; i<((D->NumarDetectii) - (Template->NumarDetectii)); i++)
    {
        (D->V+i)->valoare_corelatie=(aux+i)->valoare_corelatie;
        (D->V+i)->X1=(aux+i)->X1;
        (D->V+i)->Y1=(aux+i)->Y1;
        (D->V+i)->Y2=(aux+i)->Y2;
        (D->V+i)->X2=(aux+i)->X2;
        (D->V+i)->cifra_det=(aux+i)->cifra_det;
    }
    free(aux);
    for(i=(D->NumarDetectii) - (Template->NumarDetectii); i<D->NumarDetectii; i++)
    {
        (D->V+i)->valoare_corelatie=(Template->V+i-(D->NumarDetectii - Template->NumarDetectii))->valoare_corelatie;
        (D->V+i)->X1=(Template->V+i-(D->NumarDetectii - Template->NumarDetectii))->X1;
        (D->V+i)->Y1=(Template->V+i-(D->NumarDetectii - Template->NumarDetectii))->Y1;
        (D->V+i)->Y2=(Template->V+i-(D->NumarDetectii - Template->NumarDetectii))->Y2;
        (D->V+i)->X2=(Template->V+i-(D->NumarDetectii - Template->NumarDetectii))->X2;
        (D->V+i)->cifra_det=cifra;
    }
    free(Template->V);
}

int cmp(const void *a,const void *b)
{
    double rez=((fereastra *)b)->valoare_corelatie - ((fereastra *)a)->valoare_corelatie;
    if(rez<0)
        return -1;
    else
    {
        if(rez==0)
            return 0;
        else
            return 1;
    }
}

double suprapunere(fereastra *a,fereastra *b)
{
    int i,j,A_intersectie=0,A_a,A_b;
    for(i=a->Y1; i<=a->Y2; i++)
        for(j=a->X1; j<=a->X2; j++)
            if(i>=b->Y1 && i<=b->Y2 && j>=b->X1 && j<=b->X2)
                A_intersectie++;
    A_a=(mod((a->Y2) - (a->Y1)) + 1)*(mod((a->X2) - (a->X1)) + 1);
    A_b=(mod((b->Y2) - (b->Y1)) + 1)*(mod((b->X2) - (b->X1)) + 1);
    double rez;
    rez=((double)A_intersectie)/(double)(A_a  + A_b - A_intersectie);
    return rez;
}

void desenare(char *cale_imagine,fereastra A,RGB c)
{
    unsigned int Height,Width;
    FILE *f;
    f=fopen(cale_imagine,"rb");
    if(f==NULL)
    {
        printf("Eroare la deschiderea imaginii in care trebuie executata operatia de a desena ferestrele detectate la template matching.\n");
        return;
    }
    fseek(f,18,SEEK_SET);
    fread(&Width,sizeof(unsigned int),1,f);
    fread(&Height,sizeof(unsigned int),1,f);
    fclose(f);
    RGB *V;
    V=IncarcareFormaLiniarizata(cale_imagine);
    if(V==NULL)
    {
        printf("Eroare la liniarizarea imaginii pentru a desena ferestrelele detectate la template matching.\n");
        return ;
    }
    int i,j;
    for(i=A.Y1; i<=A.Y2; i++)
    {
        j=A.X1;
        if(i>=0 && i<Height && j>=0 && j<Width)
        {
            (V+i*Width+j)->R=c.R;
            (V+i*Width+j)->G=c.G;
            (V+i*Width+j)->B=c.B;
        }
        j=A.X2;
        if(i>=0 && i<Height && j>=0 && j<Width)
        {
            (V+i*Width+j)->R=c.R;
            (V+i*Width+j)->G=c.G;
            (V+i*Width+j)->B=c.B;
        }
    }
    for(j=A.X1; j<=A.X2; j++)
    {
        i=A.Y1;
        if(i>=0 && i<Height && j>=0 && j<Width)
        {
            (V+i*Width+j)->R=c.R;
            (V+i*Width+j)->G=c.G;
            (V+i*Width+j)->B=c.B;
        }
        i=A.Y2;
        if(i>=0 && i<Height && j>=0 && j<Width)
        {
            (V+i*Width+j)->R=c.R;
            (V+i*Width+j)->G=c.G;
            (V+i*Width+j)->B=c.B;
        }
    }
    SalvareExternFormaLiniarizata(cale_imagine,cale_imagine,V);
    free(V);
}

void eliminare_nonmaxime(ferestre_sablon *D)
{
    double prag=0.2;
    int i,j;
    qsort(D->V,D->NumarDetectii,sizeof(fereastra),cmp);
    for(i=0; i<(D->NumarDetectii) - 1; i++)
        for(j=i+1; j<D->NumarDetectii; j++)
            if((D->V+i)->valoare_corelatie!=-2 && suprapunere((D->V+i),(D->V+j))>prag)
                (D->V+j)->valoare_corelatie=-2;
}

int main()
{
    char *cale_imagine,*cale_imagine_destinatie,*cale_cheie_secreta;
    cale_imagine=(char*)malloc(100*sizeof(char));
    if(cale_imagine==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice pentru calea imaginii.\n");
        return 0;
    }
    cale_imagine_destinatie=(char*)malloc(100*sizeof(char));
    if(cale_imagine_destinatie==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice pentru calea imaginii destinatie.\n");
        return 0;
    }
    cale_cheie_secreta=(char*)malloc(100*sizeof(char));
    if(cale_cheie_secreta==NULL)
    {
        printf("Eroare la alocarea memoriei dinamice pentru calea cheii secrete.\n");
        return 0;
    }
    ferestre_sablon D;
    D.NumarDetectii=0;

    int opt;
    do
    {
        printf("******************************************************************************************************\n");
        printf("1. Cripteaza o imagine color BMP si o salveaza in memorie externa.\n");
        printf("2. Decripteaza o imagine color BMP criptata si o salveaza in memorie externa.\n");
        printf("3. Afiseaza valorile testului chi-patrat pentru imaginea initiala si cea criptata.\n");
        printf("4. Ruleaza operatia de template matching petru o imagine color BMP si o colectie de sabloane color BMP.\n");
        printf("5. Ruleaza functia de eliminare a non-maximelor si deseneaza intr-o imagine copie a primei imagini pentru a o nu altera detectiile ramase cu o culoare specifica.\n");
        printf("6. Terminare program.\n");
        printf("Alegerea dumneavoastra:   ");
        scanf("%d",&opt);
        switch(opt)
        {
        case 1:
        {
            printf("\nDati calea imaginii ce doriti sa fie criptata:   ");
            scanf("%s",cale_imagine);
            printf("\nDati calea imaginii destinatie a imiaginii criptate:   ");
            scanf("%s",cale_imagine_destinatie);
            printf("\nDati calea fisierului cu cheia secreta:   ");
            scanf("%s",cale_cheie_secreta);
            criptare(cale_imagine,cale_imagine_destinatie,cale_cheie_secreta);
            break;
        }
        case 2:
        {
            printf("\nDati calea imaginii criptate ce doriti sa fie decriptata:   ");
            scanf("%s",cale_imagine);
            printf("\nDati calea imaginii destinatie a imiaginii decriptate:   ");
            scanf("%s",cale_imagine_destinatie);
            printf("\nDati calea fisierului cu cheia secreta:   ");
            scanf("%s",cale_cheie_secreta);
            decriptare(cale_imagine,cale_imagine_destinatie,cale_cheie_secreta);
            break;
        }
        case 3:
        {
            printf("\nDati calea imaginii decriptate:   ");
            scanf("%s",cale_imagine);
            test_chi_patrat(cale_imagine);
            printf("\nDati calea imaginii criptate:   ");
            scanf("%s",cale_imagine);
            test_chi_patrat(cale_imagine);
            break;
        }
        case 4:
        {
            double ps=0.49; //am pus acest prag deoarece la 0.5 nu imi detecta primul 3
            ferestre_sablon Template;
            Template.NumarDetectii=0;
            FILE *f,*g;
            f=fopen("Sabloane.txt","r");
            if(f==NULL)
            {
                printf("Eroare la deschiderea fisierului cu sabloane si imaginea pentru template matching.\n");
                return 0;
            }
            g=fopen("Sabloane_grayscale.txt","w");
            if(g==NULL)
            {
                printf("Eroare la deschiderea fisierului pentru a scrie denumrile sabloanelor grayscale.\n");
                return 0;
            }

            //transformrea imaginii si a sabloanelor in imagini grayscale la destinatiile lor initiale cu sufixul _grayscale
            while(1)
            {
                fscanf(f,"%s",cale_imagine);
                if(feof(f))
                    break;
                strncpy(cale_imagine_destinatie,cale_imagine,strlen(cale_imagine)-4);
                cale_imagine_destinatie[strlen(cale_imagine)-4]=0;
                strcat(cale_imagine_destinatie,"_grayscale.bmp");
                cale_imagine_destinatie[strlen(cale_imagine_destinatie)]=0;
                grayscale_image(cale_imagine,cale_imagine_destinatie);
                fprintf(g,"%s\n",cale_imagine_destinatie);
            }
            fclose(f);
            fclose(g);

            f=fopen("Sabloane_grayscale.txt","r");
            if(f==NULL)
            {
                printf("Eroare la deschiderea fisierului pentru a scrie denumrile sabloanelor grayscale.\n");
                return 0;
            }
            fscanf(f,"%s",cale_imagine);
            int cifra=0;
            while(1)
            {
                fscanf(f,"%s",cale_imagine_destinatie);
                if(feof(f))
                    break;

                Template=template_matching(cale_imagine,cale_imagine_destinatie,ps);
                if(Template.NumarDetectii==-1) //a avut o eroare in functia de template_matching
                    return 0;
                creare_tablou(&D,&Template,cifra);
                cifra++;
            }
            fclose(f);
            break;
        }
        case 5:
        {
            if(D.NumarDetectii)
            {
                int i;
                eliminare_nonmaxime(&D);
                RGB *Culori;
                Culori=(RGB*) malloc(10*sizeof(RGB));
                if(Culori==NULL)
                {
                    printf("Eroare la alocarea memoriei dinamice pentru vectorul de structuri al culorilor cu care se deseneaza ferestrele detectate.\n");
                    return 0;
                }
                FILE *f=fopen("culori.txt","r");
                if(f==NULL)
                {
                    printf("Eroare la deschiderea fisierului ce contine valorile RGB ale culorilor.\n");
                    return 0;
                }
                unsigned int r,g,b;
                unsigned char *aux;
                for(i=0; i<=9; i++)
                {
                    fscanf(f,"%u",&r);
                    fscanf(f,"%u",&g);
                    fscanf(f,"%u",&b);
                    aux=(unsigned char*)&r;
                    (Culori+i)->R=*aux;
                    aux=(unsigned char*)&g;
                    (Culori+i)->G=*aux;
                    aux=(unsigned char*)&b;
                    (Culori+i)->B=*aux;

                }
                fclose(f);
                f=fopen("Sabloane.txt","r");
                if(f==NULL)
                {
                    printf("Eroare la deschiderea fisierului ce contine numele sabloanelor.\n");
                    return 0;
                }
                fscanf(f,"%s",cale_imagine);
                fclose(f);
                //cream copia imaginii peste care vom desena
                strncpy(cale_imagine_destinatie,cale_imagine,strlen(cale_imagine)-4);
                cale_imagine_destinatie[strlen(cale_imagine)-4]=0;
                strcat(cale_imagine_destinatie,"_desenata.bmp");
                cale_imagine_destinatie[strlen(cale_imagine_destinatie)]=0;
                SalvareExternFormaLiniarizata(cale_imagine,cale_imagine_destinatie,IncarcareFormaLiniarizata(cale_imagine));
                for(i=0; i<D.NumarDetectii; i++)
                    if((D.V+i)->valoare_corelatie > 0)
                        desenare(cale_imagine_destinatie,*(D.V+i),*(Culori+((D.V+i)->cifra_det)));
                free(D.V);
                free(Culori);
            }
            else
            {
                printf("Nu exista detectii sau nu s-a efectuat operatia de template matching.\n");
            }
            break;
        }
        default:
            free(cale_imagine);
            free(cale_cheie_secreta);
            free(cale_imagine_destinatie);
            break;
        }
    }
    while(opt<6 && opt>=1);
    return 0;
}
