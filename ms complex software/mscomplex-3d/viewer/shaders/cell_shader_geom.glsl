#version 120
#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable

//HEADER_REPLACE_BEGIN

const float even_sz = 0.2;
const float odd_sz  = 0.6;

//HEADER_REPLACE_END

vec4 diffuse;
vec4 ambient;
vec4 in_color;

vec3 lightDir;
vec3 halfVector;

void set_front_color_xfm(vec4 p1,vec4 p2,vec4 p3)
{
  vec3 v1 = (p1.xyz-p2.xyz);
  vec3 v2 = (p1.xyz-p3.xyz);

  vec3  halfV,viewV,ldir;
  float NdotL,NdotHV;
  vec4  color  = ambient;

  vec3  n  = normalize(cross(v1,v2));
  NdotL = max(dot(n,lightDir),0.0);

  if (NdotL > 0.0)
  {
    halfV = normalize(halfVector);
    NdotHV = max(dot(n,halfV),0.0);
    if(NdotHV > 0.0)
    {
      color += gl_FrontMaterial.specular *gl_LightSource[0].specular *pow(NdotHV,gl_FrontMaterial.shininess);
      color += diffuse * NdotL;
    }
  }
  gl_FrontColor = in_color*color;
}

void set_light_constants()
{
  diffuse    = gl_LightSource[0].diffuse;
  ambient    = gl_LightSource[0].ambient + gl_LightModel.ambient;
  lightDir   = normalize(vec3(gl_LightSource[0].position));
  halfVector = normalize(gl_LightSource[0].halfVector.xyz);
}

vec3[8] get_box(vec3 c)
{
  vec3[8] p;  int pos = 0;  ivec3 i;  vec3  sz;

  sz.x = (((int(c.x))&1) == 1)?(odd_sz):(even_sz);
  sz.y = (((int(c.y))&1) == 1)?(odd_sz):(even_sz);
  sz.z = (((int(c.z))&1) == 1)?(odd_sz):(even_sz);

  for(i[2] = -1 ; i[2] <= 1 ;i[2]+=2)
    for(i[1] = -1 ; i[1] <= 1 ;i[1]+=2)
      for(i[0] = -1 ; i[0] <= 1 ;i[0]+=2)
        p[pos++] = c+i*sz;

  return p;
}

void draw_quad_post_xfm(vec4 p1,vec4 p2,vec4 p3,vec4 p4)
{
  set_front_color_xfm(p1,p2,p3);

  gl_Position = p1; EmitVertex();
  gl_Position = p2; EmitVertex();
  gl_Position = p3; EmitVertex();
  EndPrimitive();

  gl_Position = p3; EmitVertex();
  gl_Position = p2; EmitVertex();
  gl_Position = p4; EmitVertex();

  EndPrimitive();
}

void draw_box(vec3[8] b)
{
  vec4[8] xb;

  for ( int i = 0 ; i <8 ; ++i)
    xb[i] = (gl_ModelViewProjectionMatrix*vec4(b[i],1.0));

  draw_quad_post_xfm(xb[1],xb[0],xb[3],xb[2]);
  draw_quad_post_xfm(xb[7],xb[6],xb[5],xb[4]);

  draw_quad_post_xfm(xb[0],xb[1],xb[4],xb[5]);
  draw_quad_post_xfm(xb[3],xb[2],xb[7],xb[6]);

  draw_quad_post_xfm(xb[0],xb[4],xb[2],xb[6]);
  draw_quad_post_xfm(xb[5],xb[1],xb[7],xb[3]);
}

void main()
{
  set_light_constants();

  for(int i=0; i< gl_VerticesIn; i++)
  {
    in_color = gl_FrontColorIn[i];

    draw_box(get_box(gl_PositionIn[i].xyz));
  } 
}