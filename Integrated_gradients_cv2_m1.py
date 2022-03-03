
#Import libraries
import tensorflow as tf
import tensorflow.keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten, Conv3D, MaxPooling3D,BatchNormalization, SpatialDropout3D
import numpy as np

#GPU settings 
config=tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True
mirrored_strategy = tf.distribute.MirroredStrategy()

#input
input_shape=(128,128,128,1)
#make empty model
def just_make():
    with mirrored_strategy.scope():
        model=Sequential()
        model.add(Conv3D(32,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=input_shape))
        model.add(BatchNormalization())
        model.add(Conv3D(32,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=input_shape))
     #   model.add(BatchNormalization())
        model.add(SpatialDropout3D(0.2))
        model.add(MaxPooling3D(pool_size=(2,2,2)))
        model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(64,64,64)))
        model.add(BatchNormalization())
        model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(64,64,64)))
     #   model.add(BatchNormalization())
        model.add(SpatialDropout3D(0.2))
        model.add(MaxPooling3D(pool_size=(2,2,2)))
        model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(32,32,32)))
        model.add(BatchNormalization())
        model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(32,32,32)))
        model.add(BatchNormalization())
    #     model.add(SpatialDropout3D(0.2))
        model.add(MaxPooling3D(pool_size=(2,2,2)))
        model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(16,16,16)))
        model.add(BatchNormalization())
        model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(16,16,16)))
        model.add(BatchNormalization())
        model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(16,16,16)))
        model.add(BatchNormalization())
    #     model.add(SpatialDropout3D(0.2))
        model.add(MaxPooling3D(pool_size=(2,2,2)))
        model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(8,8,8)))
        model.add(BatchNormalization())
        model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(8,8,8)))
        model.add(BatchNormalization())
        model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(8,8,8)))
        model.add(BatchNormalization())
    #     model.add(SpatialDropout3D(0.2))
        model.add(MaxPooling3D(pool_size=(2,2,2)))
        
        model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(4,4,4)))
        model.add(BatchNormalization())
        model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(4,4,4)))
        model.add(BatchNormalization())
        model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(4,4,4)))
        model.add(BatchNormalization())
        
        model.add(MaxPooling3D(pool_size=(2,2,2)))
        
        model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(2,2,2)))
        model.add(BatchNormalization())
        model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(2,2,2)))
        model.add(BatchNormalization())
        model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(2,2,2)))
        model.add(BatchNormalization())
        model.add(Flatten())
        model.add(Dense(96, activation='relu', kernel_initializer='he_uniform'))
        #     model.add(Dropout(0.2))
        model.add(Dense(32, activation='relu',kernel_initializer='he_uniform'))
        model.add(Dense(1))
    return model


ens=just_make()
ens.load_weights('FINAL_oneset_CV2_m1')

def interpolate_images(baseline,image,alphas):
    alphas_x=alphas[:,tf.newaxis,tf.newaxis,tf.newaxis,tf.newaxis]
    baseline_x=tf.expand_dims(baseline,axis=0)
    baseline_x=tf.cast(baseline_x,tf.float32)

    input_x=tf.expand_dims(image,axis=0)
    input_x=tf.cast(input_x,tf.float32)
    delta=input_x-baseline_x
    delta=tf.cast(delta,tf.float32)
    images=baseline_x+alphas_x*delta
    print('i i done')
    return images

def compute_gradients(images):
    with tf.GradientTape() as tape:
        tape.watch(images)
        predicted_age=ens(images)

        #probs=tf.nn.softmax(logits,axis=-1)[:,target_class_idx]
    #print(predicted_age)
    print('cg done')
    return tape.gradient(predicted_age,images)

def integral_approximation(gradients):
    grads=(gradients[:-1]+gradients[1:])/tf.constant(2.0)
    integrated_gradients=tf.math.reduce_mean(grads,axis=0)
    print('ia done')
    return integrated_gradients

#@tf.function
def integrated_gradients(model,baseline,image,m_steps=100,batch_size=8):
    #1 Generate alphas
    alphas=tf.linspace(start=0.0, stop=1.0,num=m_steps+1)

    #Init tensor array outside loop to collect gradients
    gradient_batches = tf.TensorArray(tf.float32,size=m_steps+1)
    k=0
    #Iterate alphas range and batch computation for spee,memory efficiency and scaling to larger m steps 
    for alpha in tf.range(0,alphas.shape[0],batch_size):
        from_=alpha
        to=tf.minimum(from_ + batch_size,alphas.shape[0])
        alpha_batch=alphas[from_:to]

        #Generate interpolated inputs between baseline and input
        interpolated_path_input_batch=interpolate_images(baseline=baseline, image=image, alphas=alpha_batch)
        print(interpolated_path_input_batch.shape)
        #batch 128 128 128 1
        #compute gradients between model outputs and interpolated outputs
        gradient_batch=compute_gradients(model=model,images=interpolated_path_input_batch)
        print(gradient_batch.shape) # None 128 128 128 1 #has to be batch 128 128 128 1 
        #write batch indiceds and gradients to extend TensorArray
        gradient_batches=gradient_batches.scatter(tf.range(from_,to),gradient_batch)
        k+=1
        print(k)
    #Stack path gradients together row-wise into single tensor 
    total_gradients=gradient_batches.stack()
    print(total_gradients.shape)
    #4. Integral approximation through averaging gradients 
    avg_gradients= integral_approximation(gradients=total_gradients)
    print(avg_gradients.shape)
    #avg_gradients=tf.cast(avg_gradients,tf.float32)
    #5. Scale integrated gradients with respect to input
    #baseline=tf.expand_dims(baseline,axis=3)


    #input_x=tf.expand_dims(image,axis=3)
    pre=image-baseline
    pre=tf.cast(pre,tf.float32) # 128 128 128 1 
    print(pre.shape)
    integrated_gradients=pre*avg_gradients
    #integrated_gradients=1
    return integrated_gradients

#integrated_gradients(ens,baseline,image)
#@tf.function
def integrated_gradients2(baseline,image,m_steps=100,batch_size=8):
    # 1. Generate alphas
    alphas = tf.linspace(start=0.0, stop=1.0, num=m_steps)

    # Accumulate gradients across batches
    integrated_gradients = 0.0

    # Batch alpha images
    ds = tf.data.Dataset.from_tensor_slices(alphas).batch(batch_size)
    for batch in ds:

        # 2. Generate interpolated images
        batch_interpolated_inputs = interpolate_images(baseline=baseline,image=image,alphas=batch)
        print(batch_interpolated_inputs.shape)
        # 3. Compute gradients between model outputs and interpolated inputs
        batch_gradients = compute_gradients(images=batch_interpolated_inputs)
        print(batch_gradients.shape)
        # 4. Average integral approximation. Summing integrated gradients across batches.
        integrated_gradients += integral_approximation(gradients=batch_gradients)
        print(integrated_gradients.shape)

    # 5. Scale integrated gradients with respect to input
    scaled_integrated_gradients = (image - baseline) * integrated_gradients                                     
    return scaled_integrated_gradients

#List of 100 individuals with the youngest predicted age 
young_list=np.loadtxt('young_list.txt')
ID128=np.load('/media/leelabsg-storage0/jangho/Aging_npy/ID_128_full.npy')
first=0
for young in young_list:
    location=int(np.where(ID128==young)[0])
    quot=(location//200)+1
    mri=np.load('/media/leelabsg-storage0/jangho/Aging_npy/final_array_128_full_'+str(quot)+'.npy')
    remnants=location%200
    image=mri[:,:,:,remnants]
    image=np.expand_dims(image,axis=3)
    #image.astype(np.float16)
    baseline=np.zeros(shape=(128,128,128,1))

    result=integrated_gradients2(baseline=baseline,image=image)
    print(result.shape)
    result=result.numpy()

    if first==0:
        final=result
        first=1
    else:
        final=final+result
    print(young)
    
np.save('young_IG100_result_cv2_m1',final/100)


