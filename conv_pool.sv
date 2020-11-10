parameter BIT_WIDTH=8;
module conv_pool (
    input logic clk, rst,
    input logic [127:0] image_4x4_r,
    input logic [127:0] image_4x4_g,
    input logic [127:0] image_4x4_b,

    input logic [71:0]  kernel_r,
    input logic [71:0]  kernel_g,
    input logic [71:0]  kernel_b,
    
    output logic        input_re,
    output logic [15:0] input_addr,

    output logic        output_we,
    output logic [15:0] output_addr,

    output logic [7:0]  y
);

//////////REGISTERS/////////
//
// FLOPPED I/O
//
////////////////////////////
logic [127:0] f_image_4x4_r;
logic [127:0] f_image_4x4_g;
logic [127:0] f_image_4x4_b;
logic [71:0]  f_kernel_r; 
logic [71:0]  f_kernel_g;
logic [71:0]  f_kernel_b;
logic         f_input_re;
logic [15:0]  f_input_addr;
logic         f_output_we;
logic [15:0]  f_output_addr;

always_ff @(posedge clk) begin
	if (rst) begin
	///////////INPUTS////////////////
	    f_image_4x4_r <= image_4x4_r;   
	    f_image_4x4_g <= image_4x4_g; 
	    f_image_4x4_b <= image_4x4_b;
	    f_kernel_r    <= kernel_r;  
	    f_kernel_g    <= kernel_g;
	    f_kernel_b    <= kernel_b;
	//////////OUTPUTS///////////////
	    input_re      <= f_input_re;
	    output_we     <= f_output_we;
	    input_addr    <= f_input_addr;
	    output_addr   <= f_output_addr;
	end
	else begin end
end
///////////////////////////////////
//
// RGB IMAGE SIGN CAST
// s_f_* (SIGNED/FLOPPED)
//
////////////////////////////////////////// 
logic signed [8:0] s_f_image_4x4_r [15:0];    
logic signed [8:0] s_f_image_4x4_g [15:0];
logic signed [8:0] s_f_image_4x4_b [15:0];
logic signed [7:0] s_f_kernel_r [8:0];
logic signed [7:0] s_f_kernel_g [8:0];
logic signed [7:0] s_f_kernel_b [8:0];

genvar i,ii;
generate
    for (i = 1; i < 17; i = i+1) begin
        assign s_f_image_4x4_r[i-1] = f_image_4x4_r[(i*BIT_WIDTH)-1:(i*BIT_WIDTH)-8];
        assign s_f_image_4x4_g[i-1] = f_image_4x4_g[(i*BIT_WIDTH)-1:(i*BIT_WIDTH)-8];
        assign s_f_image_4x4_b[i-1] = f_image_4x4_b[(i*BIT_WIDTH)-1:(i*BIT_WIDTH)-8];
    end

    for (ii = 1; ii < 10; ii = ii+1) begin
        assign s_f_kernel_r[ii-1] = f_kernel_r[(ii*BIT_WIDTH)-1:(ii*BIT_WIDTH)-8];
        assign s_f_kernel_g[ii-1] = f_kernel_g[(ii*BIT_WIDTH)-1:(ii*BIT_WIDTH)-8];
        assign s_f_kernel_b[ii-1] = f_kernel_b[(ii*BIT_WIDTH)-1:(ii*BIT_WIDTH)-8];
    end
endgenerate 
//////////////////////////////////////////////////////////////////////////////////
// INTERMEDIATE MULTS FOR 3x3 KERNEL WITH 4x4 IMAGE
// mult_ITER_COLOR [INDX] ***SIGNED 17'b 
// 4 ITERATIONS for each OP
// 9 INDICIES since 3x3 mask
// 3 COLORS
//
///////////////////// RED ////////////////////
logic signed [(2*BIT_WIDTH):0] mult_1_r [8:0];
logic signed [(2*BIT_WIDTH):0] mult_2_r [8:0];
logic signed [(2*BIT_WIDTH):0] mult_3_r [8:0];
logic signed [(2*BIT_WIDTH):0] mult_4_r [8:0];
/////////////////// GREEN ////////////////////
logic signed [(2*BIT_WIDTH):0] mult_1_g [8:0];
logic signed [(2*BIT_WIDTH):0] mult_2_g [8:0];
logic signed [(2*BIT_WIDTH):0] mult_3_g [8:0];
logic signed [(2*BIT_WIDTH):0] mult_4_g [8:0];
/////////////////// BLUE /////////////////////
logic signed [(2*BIT_WIDTH):0] mult_1_b [8:0];
logic signed [(2*BIT_WIDTH):0] mult_2_b [8:0];
logic signed [(2*BIT_WIDTH):0] mult_3_b [8:0];
logic signed [(2*BIT_WIDTH):0] mult_4_b [8:0];
//////////////////////////////////////////////
//
// COMB LOGIC
//
//////////////////////////////////////////////
//
genvar j;
generate
    for (j = 0; j < 3; j = j+1) begin
        /////////////////////////// RED ///////////////////////////////
        assign mult_1_r[j]   = s_f_image_4x4_r[j]    * s_f_kernel_r[j]; 
        assign mult_1_r[j+3] = s_f_image_4x4_r[j+4]  * s_f_kernel_r[j+3]; 
        assign mult_1_r[j+6] = s_f_image_4x4_r[j+8]  * s_f_kernel_r[j+6]; 
         
        assign mult_2_r[j]   = s_f_image_4x4_r[j+1]  * s_f_kernel_r[j]; 
        assign mult_2_r[j+3] = s_f_image_4x4_r[j+5]  * s_f_kernel_r[j+3]; 
        assign mult_2_r[j+6] = s_f_image_4x4_r[j+9]  * s_f_kernel_r[j+6];
       
        assign mult_3_r[j]   = s_f_image_4x4_r[j+4]  * s_f_kernel_r[j]; 
        assign mult_3_r[j+3] = s_f_image_4x4_r[j+8]  * s_f_kernel_r[j+3]; 
        assign mult_3_r[j+6] = s_f_image_4x4_r[j+12] * s_f_kernel_r[j+6];

        assign mult_4_r[j]   = s_f_image_4x4_r[j+5]  * s_f_kernel_r[j]; 
        assign mult_4_r[j+3] = s_f_image_4x4_r[j+9]  * s_f_kernel_r[j+3]; 
        assign mult_4_r[j+6] = s_f_image_4x4_r[j+13] * s_f_kernel_r[j+6];
        ///////////////////////////// GREEN ////j////////////////////////
        assign mult_1_g[j]   = s_f_image_4x4_g[j]    * s_f_kernel_g[j]; 
        assign mult_1_g[j+3] = s_f_image_4x4_g[j+4]  * s_f_kernel_g[j+3]; 
        assign mult_1_g[j+6] = s_f_image_4x4_g[j+8]  * s_f_kernel_g[j+6]; 
         
        assign mult_2_g[j]   = s_f_image_4x4_g[j+1]  * s_f_kernel_g[j]; 
        assign mult_2_g[j+3] = s_f_image_4x4_g[j+5]  * s_f_kernel_g[j+3]; 
        assign mult_2_g[j+6] = s_f_image_4x4_g[j+9]  * s_f_kernel_g[j+6];
       
        assign mult_3_g[j]   = s_f_image_4x4_g[j+4]  * s_f_kernel_g[j]; 
        assign mult_3_g[j+3] = s_f_image_4x4_g[j+8]  * s_f_kernel_g[j+3]; 
        assign mult_3_g[j+6] = s_f_image_4x4_g[j+12] * s_f_kernel_g[j+6];

        assign mult_4_g[j]   = s_f_image_4x4_g[j+5]  * s_f_kernel_g[j]; 
        assign mult_4_g[j+3] = s_f_image_4x4_g[j+9]  * s_f_kernel_g[j+3]; 
        assign mult_4_g[j+6] = s_f_image_4x4_g[j+13] * s_f_kernel_g[j+6];
        //////////////////////////// BLUE /////j/////////////////////////
        assign mult_1_b[j]   = s_f_image_4x4_b[j]    * s_f_kernel_b[j]; 
        assign mult_1_b[j+3] = s_f_image_4x4_b[j+4]  * s_f_kernel_b[j+3]; 
        assign mult_1_b[j+6] = s_f_image_4x4_b[j+8]  * s_f_kernel_b[j+6]; 
         
        assign mult_2_b[j]   = s_f_image_4x4_b[j+1]  * s_f_kernel_b[j]; 
        assign mult_2_b[j+3] = s_f_image_4x4_b[j+5]  * s_f_kernel_b[j+3]; 
        assign mult_2_b[j+6] = s_f_image_4x4_b[j+9]  * s_f_kernel_b[j+6];
       
        assign mult_3_b[j]   = s_f_image_4x4_b[j+4]  * s_f_kernel_b[j]; 
        assign mult_3_b[j+3] = s_f_image_4x4_b[j+8]  * s_f_kernel_b[j+3]; 
        assign mult_3_b[j+6] = s_f_image_4x4_b[j+12] * s_f_kernel_b[j+6];

        assign mult_4_b[j]   = s_f_image_4x4_b[j+5]  * s_f_kernel_b[j]; 
        assign mult_4_b[j+3] = s_f_image_4x4_b[j+9]  * s_f_kernel_b[j+3]; 
        assign mult_4_b[j+6] = s_f_image_4x4_b[j+13] * s_f_kernel_b[j+6];
end
endgenerate
//////////////////////////////////////////////
//
// 2x2x3 ACCUMULATE REGISTERS
// acc_COLOR_ITER ***SIGNED 21'b (+ 4 BITS FOR OVRFLW)
// 3 COLORS
// 4 ACCUMULATION ITERATIONS ---> 2x2x3
//
////////////////////////////////
logic signed [20:0] acc_r [3:0]; 
logic signed [20:0] acc_g [3:0]; 
logic signed [20:0] acc_b [3:0]; 

always @ (posedge clk) begin
acc_r[0] <= (mult_1_r[0] + mult_1_r[1] + mult_1_r[2] + mult_1_r[3] + mult_1_r[4] + mult_1_r[5] + mult_1_r[6] + mult_1_r[7] + mult_1_r[8]);                            
acc_r[1] <= (mult_2_r[0] + mult_2_r[1] + mult_2_r[2] + mult_2_r[3] + mult_2_r[4] + mult_2_r[5] + mult_2_r[6] + mult_2_r[7] + mult_2_r[8]);                            
acc_r[2] <= (mult_3_r[0] + mult_3_r[1] + mult_3_r[2] + mult_3_r[3] + mult_3_r[4] + mult_3_r[5] + mult_3_r[6] + mult_3_r[7] + mult_3_r[8]);                            
acc_r[3] <= (mult_4_r[0] + mult_4_r[1] + mult_4_r[2] + mult_4_r[3] + mult_4_r[4] + mult_4_r[5] + mult_4_r[6] + mult_4_r[7] + mult_4_r[8]);                            

acc_g[0] <= (mult_1_g[0] + mult_1_g[1] + mult_1_g[2] + mult_1_g[3] + mult_1_g[4] + mult_1_g[5] + mult_1_g[6] + mult_1_g[7] + mult_1_g[8]);                            
acc_g[1] <= (mult_2_g[0] + mult_2_g[1] + mult_2_g[2] + mult_2_g[3] + mult_2_g[4] + mult_2_g[5] + mult_2_g[6] + mult_2_g[7] + mult_2_g[8]);                            
acc_g[2] <= (mult_3_g[0] + mult_3_g[1] + mult_3_g[2] + mult_3_g[3] + mult_3_g[4] + mult_3_g[5] + mult_3_g[6] + mult_3_g[7] + mult_3_g[8]);                            
acc_g[3] <= (mult_4_g[0] + mult_4_g[1] + mult_4_g[2] + mult_4_g[3] + mult_4_g[4] + mult_4_g[5] + mult_4_g[6] + mult_4_g[7] + mult_4_g[8]); 

acc_b[0] <= (mult_1_b[0] + mult_1_b[1] + mult_1_b[2] + mult_1_b[3] + mult_1_b[4] + mult_1_b[5] + mult_1_b[6] + mult_1_b[7] + mult_1_b[8]);                            
acc_b[1] <= (mult_2_b[0] + mult_2_b[1] + mult_2_b[2] + mult_2_b[3] + mult_2_b[4] + mult_2_b[5] + mult_2_b[6] + mult_2_b[7] + mult_2_b[8]);                            
acc_b[2] <= (mult_3_b[0] + mult_3_b[1] + mult_3_b[2] + mult_3_b[3] + mult_3_b[4] + mult_3_b[5] + mult_3_b[6] + mult_3_b[7] + mult_3_b[8]);                            
acc_b[3] <= (mult_4_b[0] + mult_4_b[1] + mult_4_b[2] + mult_4_b[3] + mult_4_b[4] + mult_4_b[5] + mult_4_b[6] + mult_4_b[7] + mult_4_b[8]);
end


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// 2x2x1 WIRE
//                                 
// 2x2x3 --> 2x2x1 ***UNSIGNED     
// w_2x2x1_INDX ; 4 INDICIES 
//
// &
//
// SATURATION LOGIC
//
//////////////////////////////
logic signed [22:0] w_2x2x1_0;
logic signed [22:0] w_2x2x1_1;
logic signed [22:0] w_2x2x1_2;
logic signed [22:0] w_2x2x1_3;

logic [(BIT_WIDTH-1):0] w_2x2x1_sat_0;
logic [(BIT_WIDTH-1):0] w_2x2x1_sat_1;
logic [(BIT_WIDTH-1):0] w_2x2x1_sat_2;
logic [(BIT_WIDTH-1):0] w_2x2x1_sat_3;

assign w_2x2x1_0 = acc_r[0] + acc_g[0] + acc_b[0];
assign w_2x2x1_1 = acc_r[1] + acc_g[1] + acc_b[1];
assign w_2x2x1_2 = acc_r[2] + acc_g[2] + acc_b[2];
assign w_2x2x1_3 = acc_r[3] + acc_g[3] + acc_b[3];

assign w_2x2x1_sat_0 = w_2x2x1_0 > 255 ? 255 : (w_2x2x1_0 < 0 ? 0 : w_2x2x1_0);
assign w_2x2x1_sat_1 = w_2x2x1_1 > 255 ? 255 : (w_2x2x1_1 < 0 ? 0 : w_2x2x1_1);
assign w_2x2x1_sat_2 = w_2x2x1_2 > 255 ? 255 : (w_2x2x1_2 < 0 ? 0 : w_2x2x1_2);
assign w_2x2x1_sat_3 = w_2x2x1_3 > 255 ? 255 : (w_2x2x1_3 < 0 ? 0 : w_2x2x1_3); 
///////////////////////////////////////////////////////////////////////////////
//
// OUTPUT LOGIC
//
//////////////////////////////
logic [9:0] y_nxt_tmp;
logic unsigned [15:0] cnt;
logic unsigned [15:0] w_cnt;

assign y_nxt_tmp = (w_2x2x1_sat_0 + w_2x2x1_sat_1 + w_2x2x1_sat_2 + w_2x2x1_sat_3) >> 2;

logic [1:0] state,state_nxt;

parameter ST0 = 0;
parameter ST1 = 1;
parameter ST2 = 2;

always @ (posedge clk or negedge rst)
begin
	if (~rst) state <= ST0;
	else	  state <= state_nxt;
end

assign cnt_sat = cnt > 65024 ? 1 : 0;

always @ (posedge clk)
begin
	if (~rst) begin
		state_nxt     <= ST0;
		cnt 	      <= 16'b0;
		f_input_re    <=  1'b0;
		f_input_addr  <= 16'b0;
		
	end
	else begin
	case (state)
            ST0: begin 
		if (rst && cnt < 65030) begin
			if      (!cnt_sat) begin f_input_re <= 1'b1; end
			else if ( cnt_sat) begin f_input_re <= 1'b0; end

			f_input_addr <= cnt;
			cnt  	     <= cnt + 1'd1;
	        state_nxt    <= ST0;
		end
    end
	default: begin end
    endcase 
	end
end

always @ (posedge clk) begin
	if(~rst) begin
		f_output_we   <=  1'b0;
		f_output_addr <= 16'b0;
		w_cnt	      <= 16'b0;
	end

    if(f_input_addr > 2 && cnt < 65030) begin 
        f_output_we   <= 1'b1;
        f_output_addr <= w_cnt;
        w_cnt         <= w_cnt + 1'b1; 
    end
    else f_output_we <= 1'b0;	

    y <= y_nxt_tmp[7:0];
end

endmodule
